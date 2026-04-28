#!.venv/bin/python

import argparse
import os
import sys

import pandas as pd


def build_exist_dataframe(exist_rows):
    return pd.DataFrame(
        exist_rows,
        columns=["isoform", "junction_num", "genomic_start_coord", "genomic_end_coord"],
    )


def gtex_lr_exist(
    lr_gtex_junction_gtf, tumor_specific_junctions, junctions_matrix, output_path
):
    gtex_lr_rows = []
    for _, tumor_specific_junction in tumor_specific_junctions.iterrows():
        sj = junctions_matrix[
            (junctions_matrix["isoform"] == tumor_specific_junction["isoform"])
            & (
                junctions_matrix["junction_number"]
                == tumor_specific_junction["junction_num"]
            )
        ]
        if sj.empty:
            continue

        tmp = lr_gtex_junction_gtf[
            (lr_gtex_junction_gtf["Chromosome"] == sj["chrom"].values[0])
            & (lr_gtex_junction_gtf["Start"] == sj["genomic_start_coord"].values[0])
            & (lr_gtex_junction_gtf["End"] == sj["genomic_end_coord"].values[0])
        ]
        if not tmp.empty:
            gtex_lr_rows.append(
                [
                    tumor_specific_junction["isoform"],
                    tumor_specific_junction["junction_num"],
                    sj["genomic_start_coord"].values[0],
                    sj["genomic_end_coord"].values[0],
                ]
            )

    gtex_lr_list = build_exist_dataframe(gtex_lr_rows)
    gtex_lr_list.to_csv(output_path, sep="\t", index=False)
    return gtex_lr_list


def gtex_exist(
    gtex_junction_matrix, tumor_specific_junctions, junctions_matrix, output_path
):
    gtex_rows = []
    for _, tumor_specific_junction in tumor_specific_junctions.iterrows():
        sj = junctions_matrix[
            (junctions_matrix["isoform"] == tumor_specific_junction["isoform"])
            & (
                junctions_matrix["junction_number"]
                == tumor_specific_junction["junction_num"]
            )
        ]
        if sj.empty:
            continue

        tmp = gtex_junction_matrix[
            (gtex_junction_matrix["Chromosome"] == sj["chrom"].values[0])
            & (gtex_junction_matrix["Genome_Start"] == sj["genomic_start_coord"].values[0])
            & (gtex_junction_matrix["Genome_End"] == sj["genomic_end_coord"].values[0])
        ]
        if not tmp.empty:
            gtex_rows.append(
                [
                    tumor_specific_junction["isoform"],
                    tumor_specific_junction["junction_num"],
                    sj["genomic_start_coord"].values[0],
                    sj["genomic_end_coord"].values[0],
                ]
            )

    gtex_list = build_exist_dataframe(gtex_rows)
    gtex_list.to_csv(output_path, sep="\t", index=False)
    return gtex_list


def filter_tumor_specific_junctions(
    tumor_specific_junctions, gtex_sr_hits, gtex_lr_hits, output_path
):
    gtex_hits = set()

    if not gtex_sr_hits.empty:
        gtex_hits.update(
            tuple(row)
            for row in gtex_sr_hits[["isoform", "junction_num"]]
            .drop_duplicates()
            .itertuples(index=False, name=None)
        )
    if not gtex_lr_hits.empty:
        gtex_hits.update(
            tuple(row)
            for row in gtex_lr_hits[["isoform", "junction_num"]]
            .drop_duplicates()
            .itertuples(index=False, name=None)
        )

    if not gtex_hits:
        filtered = tumor_specific_junctions.copy()
    else:
        junction_ids = tumor_specific_junctions[["isoform", "junction_num"]].apply(
            lambda row: (row["isoform"], row["junction_num"]), axis=1
        )
        filtered = tumor_specific_junctions.loc[~junction_ids.isin(gtex_hits)].copy()

    filtered.to_csv(output_path, sep="\t", index=False)
    return filtered


def derive_filtered_output_path(tumor_specific_junctions, output_dir):
    file_name = os.path.basename(tumor_specific_junctions)
    if file_name.lower().endswith(".csv"):
        base_name = file_name[:-4]
    else:
        base_name = file_name
    return os.path.join(output_dir, f"{base_name}_gtex_filter.csv")


def validate_inputs(parser, args):
    if not os.path.isfile(args.gtex_sr_file):
        parser.error(f"GTEx short-read file not found: {args.gtex_sr_file}")
    if not os.path.isfile(args.gtex_lr_file):
        parser.error(f"GTEx long-read file not found: {args.gtex_lr_file}")
    if not os.path.isfile(args.junction_file):
        parser.error(f"Junction file not found: {args.junction_file}")
    if not os.path.isfile(args.tumor_specific_junctions):
        parser.error(
            "Tumor-specific junction file not found: "
            + args.tumor_specific_junctions
        )


def load_gtex_short_read_junctions(gtex_sr_file):
    gtex_junction_matrix = pd.read_csv(
        gtex_sr_file, skiprows=2, sep="\t", usecols=["Name"]
    )
    gtex_junction_matrix[["Chromosome", "Genome_Start", "Genome_End"]] = (
        gtex_junction_matrix["Name"].str.extract(r"(chr[^_]+)_([0-9]+)_([0-9]+)")
    )
    gtex_junction_matrix = gtex_junction_matrix.dropna(
        subset=["Chromosome", "Genome_Start", "Genome_End"]
    ).copy()
    gtex_junction_matrix["Genome_Start"] = gtex_junction_matrix["Genome_Start"].astype(
        int
    )
    gtex_junction_matrix["Genome_End"] = gtex_junction_matrix["Genome_End"].astype(int)
    return gtex_junction_matrix


def load_gtex_long_read_junctions(gtex_lr_file):
    lr_gtex_junction_gtf = pd.read_csv(
        gtex_lr_file,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 3, 4],
        names=["Chromosome", "Start", "End"],
    )
    lr_gtex_junction_gtf["Start"] = lr_gtex_junction_gtf["Start"].astype(int) - 1
    lr_gtex_junction_gtf["End"] = lr_gtex_junction_gtf["End"].astype(int)
    return lr_gtex_junction_gtf


def clean_cli_args(argv):
    return [arg for arg in argv if not str(arg).lstrip().startswith("#")]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter tumor-specific junctions absent from GTEx cohort"
    )
    parser.add_argument(
        "--gtex_sr_file",
        default="example/GTEx/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct",
        help="Path to GTEx short-read junction GCT file",
    )
    parser.add_argument(
        "--gtex_lr_file",
        default="example/GTEx/flair_filter_transcripts.gtf",
        help="Path to GTEx long-read junction GTF file",
    )
    parser.add_argument(
        "--junction_file",
        default="example/example_junctions.txt",
        help="Path to SQANTI3 junction file",
    )
    parser.add_argument(
        "--tumor_specific_junctions",
        default="output/tumor_specific_novel_junctions.csv",
        help="Path to tumor-specific novel junctions CSV from step 1",
    )
    parser.add_argument(
        "--output_dir",
        default="./output",
        help="Output directory for results",
    )
    args = parser.parse_args(clean_cli_args(sys.argv[1:]))

    validate_inputs(parser, args)
    os.makedirs(args.output_dir, exist_ok=True)

    gtex_junction_matrix = load_gtex_short_read_junctions(args.gtex_sr_file)
    lr_gtex_junction_gtf = load_gtex_long_read_junctions(args.gtex_lr_file)
    junctions_matrix = pd.read_csv(args.junction_file, sep="\t")
    tumor_specific_junctions = pd.read_csv(args.tumor_specific_junctions, sep="\t")

    gtex_sr_output = os.path.join(args.output_dir, "GTEx_cohort_exist.csv")
    gtex_lr_output = os.path.join(args.output_dir, "GTEx_LR_cohort_exist.csv")
    filtered_output = derive_filtered_output_path(
        args.tumor_specific_junctions, args.output_dir
    )

    gtex_sr_hits = gtex_exist(
        gtex_junction_matrix, tumor_specific_junctions, junctions_matrix, gtex_sr_output
    )
    gtex_lr_hits = gtex_lr_exist(
        lr_gtex_junction_gtf, tumor_specific_junctions, junctions_matrix, gtex_lr_output
    )
    filter_tumor_specific_junctions(
        tumor_specific_junctions, gtex_sr_hits, gtex_lr_hits, filtered_output
    )
