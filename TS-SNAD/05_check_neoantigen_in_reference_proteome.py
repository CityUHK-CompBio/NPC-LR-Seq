#! ./.venv/bin/python

from __future__ import annotations

import argparse
import csv
import gzip
import io
import zipfile
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Sequence, Set, TextIO, Tuple


DEFAULT_NEOJUNCTIONS = Path("data/neo-junction_neoantigen.txt")
DEFAULT_REFERENCE_FASTA = Path(
    "data/uniprotkb_Human_AND_model_organism_9606_2024_07_25_sp.fasta"
)
print(f"Default neo-junctions file: {DEFAULT_NEOJUNCTIONS}")
print(f"Default reference FASTA: {DEFAULT_REFERENCE_FASTA}")

GZIP_MAGIC = b"\x1f\x8b"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Check whether peptides from a tab-delimited file occur in a "
            "reference protein FASTA."
        )
    )
    parser.add_argument(
        "neojunctions_file",
        type=Path,
        nargs="?",
        default=DEFAULT_NEOJUNCTIONS,
        help=f"Tab-delimited input file. Default: {DEFAULT_NEOJUNCTIONS}",
    )
    parser.add_argument(
        "reference_fasta",
        type=Path,
        nargs="?",
        default=DEFAULT_REFERENCE_FASTA,
        help=(
            "Reference protein FASTA file. Supports plain text FASTA, gzip, "
            f"and single-file zip archives. Default: {DEFAULT_REFERENCE_FASTA}"
        ),
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("output/neoantigen_reference_check"),
        help="Output prefix for the generated TSV reports.",
    )
    parser.add_argument(
        "--peptide-column-name",
        default=None,
        help="Peptide column name in the neojunction table. Default: new_peptide.",
    )
    parser.add_argument(
        "--peptide-column",
        type=int,
        default=None,
        help=(
            "1-based peptide column index in the neojunction table. "
            "Used only when --peptide-column-name is not provided."
        ),
    )
    return parser.parse_args()


@contextmanager
def open_text(path: Path) -> Iterator[TextIO]:
    with path.open("rb") as probe_handle:
        magic = probe_handle.read(4)

    if magic.startswith(GZIP_MAGIC):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            yield handle
        return

    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            members = [member for member in archive.infolist() if not member.is_dir()]
            if not members:
                raise ValueError(f"Zip archive is empty: {path}")
            if len(members) > 1:
                member_names = ", ".join(member.filename for member in members)
                raise ValueError(
                    f"Zip archive must contain exactly one file: {path} ({member_names})"
                )
            with archive.open(members[0], "r") as raw_handle:
                with io.TextIOWrapper(
                    raw_handle, encoding="utf-8", errors="replace"
                ) as handle:
                    yield handle
        return

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        yield handle


def resolve_peptide_column_index(
    header: Sequence[str],
    peptide_column_name: str | None,
    peptide_column: int | None,
) -> int:
    normalized_header = [column.strip() for column in header]

    if peptide_column_name:
        try:
            return normalized_header.index(peptide_column_name.strip())
        except ValueError as exc:
            raise ValueError(
                f"Column '{peptide_column_name}' not found in header: {normalized_header}"
            ) from exc

    if peptide_column is None:
        raise ValueError("Either --peptide-column-name or --peptide-column must be set.")
    if peptide_column < 1:
        raise ValueError("--peptide-column must be >= 1")

    peptide_idx = peptide_column - 1
    if peptide_idx >= len(header):
        raise ValueError(
            f"Requested peptide column {peptide_column}, but header only has {len(header)} columns."
        )
    return peptide_idx


def load_neojunction_rows(
    path: Path,
    peptide_column_name: str | None,
    peptide_column: int | None,
) -> Tuple[List[str], List[List[str]], Counter[str], Dict[str, int], int]:
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Input file is empty: {path}") from exc

        peptide_idx = resolve_peptide_column_index(
            header, peptide_column_name, peptide_column
        )

        rows: List[List[str]] = []
        peptide_counts: Counter[str] = Counter()
        peptide_first_seen: Dict[str, int] = {}

        for row_number, row in enumerate(reader, start=2):
            if not row:
                continue
            if peptide_idx >= len(row):
                raise ValueError(
                    f"Line {row_number} has {len(row)} columns, expected at least {peptide_idx + 1}."
                )

            peptide = row[peptide_idx].strip().upper()
            rows.append(row)
            if not peptide:
                continue

            peptide_counts[peptide] += 1
            peptide_first_seen.setdefault(peptide, len(peptide_first_seen))

    return header, rows, peptide_counts, peptide_first_seen, peptide_idx


def iter_fasta_records(path: Path) -> Iterator[Tuple[str, str]]:
    with open_text(path) as handle:
        header: str | None = None
        seq_parts: List[str] = []

        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts).upper()
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)

        if header is not None:
            yield header, "".join(seq_parts).upper()


def scan_reference_proteome(
    reference_fasta: Path,
    peptide_counts: Counter[str],
) -> Set[str]:
    peptides_by_length: Dict[int, Set[str]] = {}
    for peptide in peptide_counts:
        peptides_by_length.setdefault(len(peptide), set()).add(peptide)

    found_peptides: Set[str] = set()
    remaining_peptides = set(peptide_counts)

    for _, sequence in iter_fasta_records(reference_fasta):
        seq_len = len(sequence)
        for kmer_length, peptide_set in peptides_by_length.items():
            if seq_len < kmer_length:
                continue
            candidate_peptides = peptide_set & remaining_peptides
            if not candidate_peptides:
                continue
            for start in range(seq_len - kmer_length + 1):
                kmer = sequence[start : start + kmer_length]
                if kmer in candidate_peptides:
                    found_peptides.add(kmer)
                    remaining_peptides.discard(kmer)
            if not remaining_peptides:
                return found_peptides

    return found_peptides


def write_unique_report(
    output_path: Path,
    peptide_counts: Counter[str],
    peptide_first_seen: Dict[str, int],
    found_peptides: Set[str],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ordered_peptides = sorted(peptide_counts, key=peptide_first_seen.__getitem__)

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "peptide",
                "length",
                "input_row_count",
                "found_in_reference",
            ]
        )
        for peptide in ordered_peptides:
            writer.writerow(
                [
                    peptide,
                    len(peptide),
                    peptide_counts[peptide],
                    "yes" if peptide in found_peptides else "no",
                ]
            )


def write_row_report(
    output_path: Path,
    header: Sequence[str],
    rows: Iterable[Sequence[str]],
    peptide_idx: int,
    found_peptides: Set[str],
    matched_only: bool = False,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(list(header) + ["found_in_reference"])

        for row in rows:
            peptide = row[peptide_idx].strip().upper() if peptide_idx < len(row) else ""
            is_matched = bool(peptide and peptide in found_peptides)
            if matched_only and not is_matched:
                continue
            writer.writerow(list(row) + ["yes" if is_matched else "no"])


def main() -> None:
    args = parse_args()

    peptide_column_name = args.peptide_column_name
    if peptide_column_name is None and args.peptide_column is None:
        peptide_column_name = "new_peptide"

    header, rows, peptide_counts, peptide_first_seen, peptide_idx = load_neojunction_rows(
        args.neojunctions_file,
        peptide_column_name,
        args.peptide_column,
    )
    found_peptides = scan_reference_proteome(args.reference_fasta, peptide_counts)

    unique_output = args.output_prefix.with_suffix(".unique.tsv")
    row_output = args.output_prefix.with_suffix(".rows.tsv")
    matched_row_output = args.output_prefix.with_suffix(".matched.rows.tsv")

    write_unique_report(
        unique_output,
        peptide_counts,
        peptide_first_seen,
        found_peptides,
    )
    write_row_report(
        row_output,
        header,
        rows,
        peptide_idx,
        found_peptides,
    )
    write_row_report(
        matched_row_output,
        header,
        rows,
        peptide_idx,
        found_peptides,
        matched_only=True,
    )

    total_rows = sum(peptide_counts.values())
    matched_unique = len(found_peptides)
    matched_rows = sum(peptide_counts[peptide] for peptide in found_peptides)

    print(f"Input rows with non-empty peptides: {total_rows}")
    print(f"Unique peptides: {len(peptide_counts)}")
    print(f"Matched unique peptides: {matched_unique}")
    print(f"Unmatched unique peptides: {len(peptide_counts) - matched_unique}")
    print(f"Matched rows: {matched_rows}")
    print(f"Unmatched rows: {total_rows - matched_rows}")
    print(f"Unique peptide report: {unique_output}")
    print(f"Row-level report: {row_output}")
    print(f"Matched row-level report: {matched_row_output}")


if __name__ == "__main__":
    main()
