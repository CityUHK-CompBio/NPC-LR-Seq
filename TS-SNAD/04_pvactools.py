#! ./.venv/bin/python

import argparse
import subprocess
from pathlib import Path

import pandas as pd


def run_cmd(cmd: list[str], log_path: Path | None = None):
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("w") as handle:
            subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, check=True)
    else:
        subprocess.run(cmd, check=True)


def normalize_hla_alleles(hla_input: str) -> str:
    path = Path(hla_input)
    if path.is_file():
        alleles = []
        for line in path.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            alleles.extend(part.strip() for part in line.split(",") if part.strip())
        source = str(path)
    else:
        alleles = [part.strip() for part in hla_input.split(",") if part.strip()]
        source = "command line"

    alleles = list(dict.fromkeys(alleles))
    if not alleles:
        raise ValueError(
            "--hla_alleles must be a comma-separated list or a txt file "
            "containing one HLA allele per line."
        )

    print(f"[hla] loaded {len(alleles)} HLA alleles from {source}")
    return ",".join(alleles)


def intersect_and_write_fasta(
    tumor_specific_junctions: str,
    neoantigen_file: str,
    fasta_path: Path,
):
    filtered_junctions = pd.read_csv(tumor_specific_junctions, sep="\t")
    for col in ("isoform", "junction_num"):
        if col not in filtered_junctions.columns:
            raise ValueError(f"Column '{col}' not found in {tumor_specific_junctions}")

    filtered_junctions["_junc_id"] = (
        filtered_junctions["junction_num"].astype(str).str.extract(r"_(\d+)$").astype(int)
    )
    keep = set(zip(filtered_junctions["isoform"], filtered_junctions["_junc_id"]))
    print(f"[intersect] step-02 surviving junctions: {len(keep)}")

    neoantigen_rows = pd.read_csv(neoantigen_file, sep="\t")
    for col in ("transcript_id", "junction_id", "new_peptide"):
        if col not in neoantigen_rows.columns:
            raise ValueError(f"Column '{col}' not found in {neoantigen_file}")

    before = len(neoantigen_rows)
    mask = neoantigen_rows.apply(
        lambda row: (row["transcript_id"], int(row["junction_id"])) in keep, axis=1
    )
    neoantigen_rows = neoantigen_rows[mask]
    print(f"[intersect] kept {len(neoantigen_rows)}/{before} 9-mer rows")

    intersected_tsv = fasta_path.parent / "neoantigen.intersected.tsv"
    intersected_tsv.parent.mkdir(parents=True, exist_ok=True)
    neoantigen_rows.to_csv(intersected_tsv, sep="\t", index=False)
    print(f"[intersect] saved intersected rows -> {intersected_tsv}")

    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    peptides = (
        neoantigen_rows["new_peptide"]
        .dropna()
        .astype(str)
        .pipe(lambda s: s[s.str.len() == 9])
        .pipe(lambda s: s[s.str.upper().apply(lambda peptide: set(peptide).issubset(valid_aa))])
        .str.upper()
        .drop_duplicates()
        .reset_index(drop=True)
    )

    if peptides.empty:
        raise ValueError("No valid 9-mer peptides found after intersection.")

    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    with fasta_path.open("w") as handle:
        for index, peptide in enumerate(peptides, start=1):
            handle.write(f">peptide_{index}\n{peptide}\n")

    print(f"[fasta] wrote {len(peptides)} unique 9-mer peptides -> {fasta_path}")
    return fasta_path


def binding_affinity(
    fasta_file: Path,
    output_directory: Path,
    hla_alleles: str,
    prefix: str,
):
    output_directory.mkdir(parents=True, exist_ok=True)
    log_file = output_directory / "pvactools.log"

    cmd = [
        "pvacbind",
        "run",
        "--blastp-path",
        "/usr/bin/blastp",
        str(fasta_file),
        prefix,
        hla_alleles,
        "NetMHCpan",
        str(output_directory),
        "-e1",
        "9",
        "-t",
        "128",
    ]
    print(f"[pvacbind] running binding affinity prediction -> {log_file}")
    run_cmd(cmd, log_path=log_file)


def filter_binding(output_directory: Path, prefix: str):
    mhc_dir = output_directory / "MHC_Class_I"
    in_tsv = mhc_dir / f"{prefix}.MHC_I.all_epitopes.tsv"
    out_tsv = mhc_dir / f"{prefix}.MHC_I.all_epitopes.filtered.tsv"

    if not in_tsv.exists():
        raise FileNotFoundError(f"pvacbind output not found: {in_tsv}")

    cmd = [
        "pvacbind",
        "binding_filter",
        "-b",
        "500",
        "-p",
        "0.5",
        "-m",
        "median",
        str(in_tsv),
        str(out_tsv),
    ]
    log_file = mhc_dir / "binding_filter.log"
    print(f"[pvacbind] filtering -> {out_tsv}")
    run_cmd(cmd, log_path=log_file)
    print(f"[done] filtered results saved to {out_tsv}")


def postprocess_filtered_epitopes(
    intersected_tsv: Path,
    filtered_epitopes_tsv: Path,
    output_tsv: Path,
):
    if not intersected_tsv.exists():
        raise FileNotFoundError(f"Intersected neoantigen TSV not found: {intersected_tsv}")
    if not filtered_epitopes_tsv.exists():
        raise FileNotFoundError(
            f"Filtered pVACbind epitopes TSV not found: {filtered_epitopes_tsv}"
        )

    epitopes = pd.read_csv(filtered_epitopes_tsv, sep="\t")
    for col in ("Epitope Seq", "HLA Allele"):
        if col not in epitopes.columns:
            raise ValueError(f"Column '{col}' not found in {filtered_epitopes_tsv}")

    intersected = pd.read_csv(intersected_tsv, sep="\t")
    required_columns = [
        "transcript_id",
        "new_peptide",
        "start_acid",
        "end_acid",
        "junction_id",
        "strand",
        "isoform_seq",
        "mrna_seq",
    ]
    for col in required_columns:
        if col not in intersected.columns:
            raise ValueError(f"Column '{col}' not found in {intersected_tsv}")

    epitope_hits = epitopes[["Epitope Seq", "HLA Allele"]].dropna().copy()
    epitope_hits["Epitope Seq"] = epitope_hits["Epitope Seq"].astype(str).str.strip()
    epitope_hits["HLA Allele"] = epitope_hits["HLA Allele"].astype(str).str.strip()
    epitope_hits = epitope_hits[
        (epitope_hits["Epitope Seq"] != "") & (epitope_hits["HLA Allele"] != "")
    ].copy()
    epitope_hits["_match_peptide"] = epitope_hits["Epitope Seq"].str.upper()
    epitope_hits = epitope_hits[["_match_peptide", "HLA Allele"]].drop_duplicates()

    intersected["_match_peptide"] = (
        intersected["new_peptide"].fillna("").astype(str).str.strip().str.upper()
    )
    matched = intersected.merge(epitope_hits, on="_match_peptide", how="inner")
    matched = matched[
        [
            "transcript_id",
            "new_peptide",
            "HLA Allele",
            "start_acid",
            "end_acid",
            "junction_id",
            "strand",
            "isoform_seq",
            "mrna_seq",
        ]
    ].copy()

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    matched.to_csv(output_tsv, sep="\t", index=False)
    print(
        "[postprocess] matched "
        f"{len(matched)}/{len(intersected)} peptide-HLA rows "
        f"from {len(epitope_hits)} filtered peptide-HLA hits -> {output_tsv}"
    )
    return output_tsv


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Predict HLA-I binding affinity for 9-mer neoantigens using pVACbind"
    )
    parser.add_argument(
        "--tumor_specific_junctions",
        required=True,
        help="Step-02 GTEx-filtered tumor-specific junction CSV (columns: isoform, junction_num)",
    )
    parser.add_argument(
        "--neoantigen_file",
        required=True,
        help="Step-03 9-mer neoantigen TSV (columns: transcript_id, junction_id, new_peptide, ...)",
    )
    parser.add_argument(
        "--hla_alleles",
        required=True,
        help=(
            "Comma-separated HLA alleles or a txt file path with one HLA allele "
            "per line, e.g. HLA-A*02:07,HLA-B*46:01,HLA-C*01:02 or "
            "example/hla_alleles.txt"
        ),
    )
    parser.add_argument(
        "--output_directory",
        required=True,
        help="Root output directory; pvacbind results go into <output_directory>/MHC_Class_I/",
    )
    parser.add_argument(
        "--prefix",
        default="NPC",
        help="Sample prefix used for pvacbind output naming (default: NPC)",
    )
    args = parser.parse_args()

    outdir = Path(args.output_directory)
    hla_alleles = normalize_hla_alleles(args.hla_alleles)

    fasta_path = outdir / "neoantigen.pep.fa"
    intersect_and_write_fasta(
        args.tumor_specific_junctions, args.neoantigen_file, fasta_path
    )
    binding_affinity(fasta_path, outdir, hla_alleles, args.prefix)
    filter_binding(outdir, args.prefix)
    postprocess_filtered_epitopes(
        outdir / "neoantigen.intersected.tsv",
        outdir / "MHC_Class_I" / f"{args.prefix}.MHC_I.all_epitopes.filtered.tsv",
        outdir / "neoantigen.intersected.filtered_epitopes.tsv",
    )
