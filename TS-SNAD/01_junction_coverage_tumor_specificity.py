#!.venv/bin/python

import argparse
import glob
import os

import pandas as pd


SJ_SUFFIX = ".SJ.out.tab"
JUNCTION_CATEGORY = "novel"
DEFAULT_TUMOR_SAMPLES = "example/tumor_samples.txt"
DEFAULT_CONTROL_SAMPLES = "example/control_samples.txt"


def normalize_sample_id(sample_name):
    sample_id = os.path.basename(str(sample_name).strip())
    if sample_id.endswith(SJ_SUFFIX):
        sample_id = sample_id[: -len(SJ_SUFFIX)]
    return sample_id


def normalize_sample_ids(sample_names):
    return [
        normalized
        for normalized in (normalize_sample_id(name) for name in sample_names)
        if normalized
    ]


def find_duplicates(values):
    seen = set()
    duplicates = set()
    for value in values:
        if value in seen:
            duplicates.add(value)
        else:
            seen.add(value)
    return sorted(duplicates)


def discover_sj_files(sj_path):
    sample_to_file = {}
    for file_path in glob.glob(os.path.join(sj_path, f"*{SJ_SUFFIX}")):
        sample_id = normalize_sample_id(file_path)
        if sample_id in sample_to_file:
            raise ValueError(
                f"Duplicate SJ.out.tab files found for sample '{sample_id}'."
            )
        sample_to_file[sample_id] = file_path
    return sample_to_file


def build_sj_lookup(sample_order, sample_to_file):
    sj_lookup = {}
    for sample_id in sample_order:
        sj_out = pd.read_csv(
            sample_to_file[sample_id], sep="\t", header=None, usecols=[0, 1, 2, 6]
        )
        sj_out = sj_out.rename(columns={0: "chrom", 1: "start", 2: "end", 6: "reads"})
        sj_lookup[sample_id] = (
            sj_out.groupby(["chrom", "start", "end"], sort=False)["reads"].sum()
        )
    return sj_lookup


def coverage(chrom, genomic_start_coord, genomic_end_coord, sample_order, sj_lookup):
    key = (chrom, genomic_start_coord, genomic_end_coord)
    return [int(sj_lookup[sample_id].get(key, 0)) for sample_id in sample_order]


def junction_coverage_matrix(junctions_matrix, sample_order, sj_lookup, output_path):
    selected_junctions = junctions_matrix[
        junctions_matrix["junction_category"] == JUNCTION_CATEGORY
    ].copy()

    rows = []
    for _, row in selected_junctions.iterrows():
        chrom = row["chrom"]
        genomic_start_coord = int(row["genomic_start_coord"])
        genomic_end_coord = int(row["genomic_end_coord"])
        coverage_values = coverage(
            chrom,
            genomic_start_coord,
            genomic_end_coord,
            sample_order,
            sj_lookup,
        )
        row_data = {
            "isoform": row["isoform"],
            "junction_num": row["junction_number"],
            "genomic_start_coord": genomic_start_coord,
            "genomic_end_coord": genomic_end_coord,
        }
        row_data.update(
            {
                sample_order[index]: coverage_values[index]
                for index in range(len(sample_order))
            }
        )
        rows.append(row_data)

    columns = [
        "isoform",
        "junction_num",
        "genomic_start_coord",
        "genomic_end_coord",
        *sample_order,
    ]
    junctions = pd.DataFrame(rows, columns=columns)
    junctions.to_csv(output_path, index=False)
    return junctions


def tumor_specific_junctions(
    junctions,
    tumor_samples,
    control_samples,
    output_path,
    fold_change=10,
    max_normal_sum=5,
):
    tumor_specific = []

    for _, junction in junctions.iterrows():
        tumor_values = pd.to_numeric(junction[tumor_samples], errors="coerce").fillna(0.0)
        control_values = pd.to_numeric(
            junction[control_samples], errors="coerce"
        ).fillna(0.0)

        mean_tumor = round(float(tumor_values.mean()), 2)
        mean_normal = round(float(control_values.mean()), 2)
        median_tumor = float(tumor_values.median())
        median_normal = float(control_values.median())
        max_tumor = float(tumor_values.max())
        max_normal = float(control_values.max())

        if (
            (mean_tumor > fold_change * mean_normal)
            and (control_values < max_normal_sum).all()
        ):
            tumor_specific.append(
                {
                    "isoform": junction["isoform"],
                    "junction_num": junction["junction_num"],
                    "mean_tumor": mean_tumor,
                    "mean_normal": mean_normal,
                    "median_tumor": median_tumor,
                    "median_normal": median_normal,
                    "max_tumor": max_tumor,
                    "max_normal": max_normal,
                }
            )

    df = pd.DataFrame(
        tumor_specific,
        columns=[
            "isoform",
            "junction_num",
            "mean_tumor",
            "mean_normal",
            "median_tumor",
            "median_normal",
            "max_tumor",
            "max_normal",
        ],
    )
    df.to_csv(output_path, sep="\t", index=False)
    return df


def parse_samples_from_txt(path, parser, arg_name):
    if not os.path.isfile(path):
        parser.error(f"{arg_name} sample file not found: {path}")
    try:
        with open(path, "r", encoding="utf-8") as handle:
            tokens = handle.read().split()
    except OSError as exc:
        parser.error(f"Failed to read {arg_name} sample file '{path}': {exc}")
    return normalize_sample_ids(tokens)


def parse_sample_tokens(tokens, parser, arg_name):
    parsed = []
    for token in tokens:
        token = str(token).strip()
        if not token:
            continue
        if token.lower().endswith(".txt"):
            if os.path.isfile(token):
                parsed.extend(parse_samples_from_txt(token, parser, arg_name))
                continue
            parser.error(f"{arg_name} sample file not found: {token}")
        parsed.extend(normalize_sample_ids([token]))
    return parsed


def validate_inputs(parser, args):
    if not os.path.isfile(args.junction_file):
        parser.error(f"Junction file not found: {args.junction_file}")
    if not os.path.isdir(args.sj_path):
        parser.error(f"SJ.out.tab directory not found: {args.sj_path}")

    try:
        sample_to_file = discover_sj_files(args.sj_path)
    except ValueError as exc:
        parser.error(str(exc))

    if not sample_to_file:
        parser.error(f"No SJ.out.tab files found in directory: {args.sj_path}")

    tumor_samples = parse_sample_tokens(args.tumor_samples, parser, "--tumor_samples")
    control_samples = parse_sample_tokens(
        args.control_samples, parser, "--control_samples"
    )

    if not tumor_samples:
        parser.error("At least one tumor sample is required.")
    if not control_samples:
        parser.error("At least one control sample is required.")

    tumor_duplicates = find_duplicates(tumor_samples)
    control_duplicates = find_duplicates(control_samples)
    if tumor_duplicates:
        parser.error(
            "Duplicate tumor sample names provided: " + ", ".join(tumor_duplicates)
        )
    if control_duplicates:
        parser.error(
            "Duplicate control sample names provided: " + ", ".join(control_duplicates)
        )

    overlap = sorted(set(tumor_samples) & set(control_samples))
    if overlap:
        parser.error(
            "Samples cannot be both tumor and control: " + ", ".join(overlap)
        )

    sample_order = tumor_samples + control_samples
    requested_samples = set(sample_order)
    discovered_samples = set(sample_to_file)

    if len(discovered_samples) != len(requested_samples):
        parser.error(
            "SJ.out.tab file count does not match requested samples: "
            f"found {len(discovered_samples)} files, "
            f"but {len(requested_samples)} samples were requested."
        )

    missing = sorted(requested_samples - discovered_samples)
    extra = sorted(discovered_samples - requested_samples)
    if missing or extra:
        details = []
        if missing:
            details.append("missing samples: " + ", ".join(missing))
        if extra:
            details.append("extra SJ files: " + ", ".join(extra))
        parser.error(
            "Sample names in requested groups and SJ.out.tab directory do not match: "
            + "; ".join(details)
        )

    return sample_order, tumor_samples, control_samples, sample_to_file


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate tumor-specific junction coverage"
    )
    parser.add_argument(
        "--junction_file",
        default="example/example_junctions.txt",
        help="Path to SQANTI3 junction file",
    )
    parser.add_argument(
        "--sj_path",
        default="example/SJ.out.tab",
        help="Path to SJ.out.tab directory",
    )
    parser.add_argument(
        "--tumor_samples",
        nargs="+",
        default=[DEFAULT_TUMOR_SAMPLES],
        help=(
            "Tumor sample names or existing .txt file paths "
            f"(default: {DEFAULT_TUMOR_SAMPLES})"
        ),
    )
    parser.add_argument(
        "--control_samples",
        nargs="+",
        default=[DEFAULT_CONTROL_SAMPLES],
        help=(
            "Control sample names or existing .txt file paths "
            f"(default: {DEFAULT_CONTROL_SAMPLES})"
        ),
    )
    parser.add_argument(
        "--output_dir",
        default="./output",
        help="Output directory for results",
    )
    parser.add_argument(
        "--fold_change",
        type=float,
        default=10,
        help="Multiplier used in mean_tumor > fold_change * mean_normal",
    )
    parser.add_argument(
        "--max_normal_sum",
        type=float,
        default=5,
        help="Max allowed normal coverage in each control sample",
    )
    args = parser.parse_args()

    sample_order, tumor_samples, control_samples, sample_to_file = validate_inputs(
        parser, args
    )
    os.makedirs(args.output_dir, exist_ok=True)

    junctions_matrix = pd.read_csv(args.junction_file, sep="\t")
    sj_lookup = build_sj_lookup(sample_order, sample_to_file)

    junction_matrix_output = os.path.join(
        args.output_dir, f"{JUNCTION_CATEGORY}_junctions_matrix.csv"
    )
    specific_output = os.path.join(
        args.output_dir, f"tumor_specific_{JUNCTION_CATEGORY}_junctions.csv"
    )

    junctions = junction_coverage_matrix(
        junctions_matrix,
        sample_order,
        sj_lookup,
        junction_matrix_output,
    )
    tumor_specific_junctions(
        junctions,
        tumor_samples,
        control_samples,
        specific_output,
        fold_change=args.fold_change,
        max_normal_sum=args.max_normal_sum,
    )
