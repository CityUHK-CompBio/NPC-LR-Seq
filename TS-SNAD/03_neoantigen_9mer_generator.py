#! ./.venv/bin/python

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pathos
from Bio import Seq
from Bio import SeqIO


def get_position_info(cds_gff):
    position_info = pd.read_csv(cds_gff, sep="\t", header=None, comment="#")
    position_info.columns = [
        "chr",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]
    position_info["transcript_id"] = position_info["attributes"].str.extract(
        r'transcript_id "(.*?)";'
    )
    position_info["gene_id"] = position_info["attributes"].str.extract(
        r'gene_id "(.*?)";'
    )
    position_info = position_info[
        ["chr", "source", "type", "start", "end", "strand", "transcript_id"]
    ]
    return position_info


def get_novel_junction(path):
    junctions = pd.read_csv(path, sep="\t")
    junctions = junctions[junctions["junction_category"] == "novel"]
    junctions["isoform"].unique()
    return junctions


def get_flags(isoform_seq, exons_in_isoform, cds_in_isoform):
    start_point = exons_in_isoform.iloc[0]["start"]
    end_point = exons_in_isoform.iloc[-1]["end"]
    length = end_point - start_point + 1
    exons_flag = np.array([0] * length)
    for i in range(len(exons_in_isoform)):
        left = exons_in_isoform.iloc[i]["start"] - start_point
        right = exons_in_isoform.iloc[i]["end"] - start_point + 1
        exons_flag[left:right] = i + 1
    cds_flag = exons_flag.copy()

    assert (cds_in_isoform.sort_values(by=["start"]) == cds_in_isoform).all().all()
    cds_left = cds_in_isoform.iloc[0]["start"] - start_point
    cds_right = cds_in_isoform.iloc[-1]["end"] - start_point + 1
    cds_flag[:cds_left] = 0
    cds_flag[cds_right:] = 0

    inner_index = np.where(exons_flag == 0)[0]
    assert np.all(cds_flag[inner_index]) == 0
    exons_flag = np.delete(exons_flag, inner_index)
    cds_flag = np.delete(cds_flag, inner_index)
    assert len(isoform_seq) == len(exons_flag)
    assert len(isoform_seq) == len(cds_flag)

    mrna_seq = isoform_seq[cds_flag.astype("bool")]
    mrna_seq = "".join(mrna_seq)
    junction_flag = cds_flag[cds_flag.astype("bool")]
    if np.all(exons_in_isoform["strand"] == "+"):
        peptide = Seq.translate(mrna_seq)
        peptide = peptide[:-1] if peptide[-1] == "*" else peptide
    elif np.all(exons_in_isoform["strand"] == "-"):
        rev_mrna_seq = Seq.reverse_complement(mrna_seq)
        peptide = Seq.translate(rev_mrna_seq)
        peptide = peptide[:-1] if peptide[-1] == "*" else peptide
        junction_flag = junction_flag[::-1]
    else:
        raise ValueError("isoform strand is not all + or -")

    return mrna_seq, peptide, exons_flag, cds_flag, junction_flag


def get_9mer_for_junc(junc_position, peptide):
    if (junc_position + 1) % 3 == 0:
        acid_pos = (junc_position + 1) // 3 - 1
        start_acid = acid_pos - 7 if acid_pos - 7 >= 0 else 0
        end_acid = acid_pos + 9 if acid_pos + 9 <= len(peptide) else len(peptide)
    else:
        acid_pos = (junc_position + 1) // 3
        start_acid = acid_pos - 8 if acid_pos - 8 >= 0 else 0
        end_acid = acid_pos + 9 if acid_pos + 9 <= len(peptide) else len(peptide)
    tmp = peptide[start_acid:end_acid]
    if len(tmp) >= 9:
        tmp = np.array(list(tmp))
        peptide_nine = np.lib.stride_tricks.sliding_window_view(tmp, 9)
        peptide_nine = ["".join(i) for i in peptide_nine]
    else:
        peptide_nine = ["slide-error"]
    return peptide_nine, tmp, start_acid, end_acid


def format_flag(flag):
    info = []
    for i in np.unique(flag):
        if i == 0:
            continue
        location1, location2 = np.where(flag == i)[0][[0, -1]]
        info.append((i, location1, location2))
    info = pd.DataFrame(data=np.array(info), columns=["id", "start", "end"])
    info.sort_values(by=["start"], inplace=True)
    result = ""
    for i in range(len(info)):
        result += f'{info.iloc[i]["id"]}:{info.iloc[i]["start"]}-{info.iloc[i]["end"]};'
    return result


def get_9mer(acids_seq, junction_flag, junction_id, strand):
    result = []
    acids_seq = np.array(list(acids_seq))
    if strand == "+":
        for junct in junction_id:
            left = np.where(junction_flag == junct)[0]
            right = np.where(junction_flag == junct + 1)[0]
            if len(left) == 0 or len(right) == 0:
                tmp = pd.DataFrame(
                    data=["junc-loc-in-start-or-end-code"], columns=["new_peptide"]
                )
            else:
                left, right = left[-1], right[0]
                assert left == right - 1
                assert junction_flag[left] == junction_flag[right] - 1
                nine_mer, start_to_end, start_acid, end_acid = get_9mer_for_junc(
                    left, acids_seq
                )
                tmp = pd.DataFrame(data=nine_mer, columns=["new_peptide"])
                tmp["start_to_end"] = ["".join(start_to_end)] * len(tmp)
                tmp["start_acid"] = start_acid + 1
                tmp["end_acid"] = end_acid
            tmp["junction_id"] = junct
            result.append(tmp)
    elif strand == "-":
        for junct in junction_id:
            left = np.where(junction_flag == junct + 1)[0]
            right = np.where(junction_flag == junct)[0]
            if len(left) == 0 or len(right) == 0:
                tmp = pd.DataFrame(
                    data=["junc-loc-in-start-or-end-code"], columns=["new_peptide"]
                )
            else:
                left, right = left[-1], right[0]
                assert left == right - 1
                assert junction_flag[left] == junction_flag[right] + 1
                nine_mer, start_to_end, start_acid, end_acid = get_9mer_for_junc(
                    left, acids_seq
                )
                tmp = pd.DataFrame(data=nine_mer, columns=["new_peptide"])
                tmp["start_to_end"] = ["".join(start_to_end)] * len(tmp)
                tmp["start_acid"] = start_acid + 1
                tmp["end_acid"] = end_acid
            tmp["junction_id"] = junct
            result.append(tmp)
    else:
        raise ValueError("strand must be + or -")
    result = pd.concat(result, axis=0, ignore_index=True)
    return result


def run_one_transcript(
    transcript_id, junction_id, isoform_seq, acids_seq, exons_info, cds_info
):
    strand = exons_info.iloc[0]["strand"]
    if strand == "-":
        isoform_seq = Seq.reverse_complement(isoform_seq)
    mrna_seq, acids, exons_flag, cds_flag, junc_flag = get_flags(
        np.array(list(isoform_seq)), exons_info, cds_info
    )
    if acids != acids_seq:
        print("transcript id", transcript_id)
        print("acids", acids)
        print("acids_seq", acids_seq)
        return None
    df_9mer = get_9mer(acids, junc_flag, junction_id, strand)
    df_9mer.insert(0, "transcript_id", transcript_id)
    df_9mer["strand"] = strand
    full_info = {
        "isoform_seq": isoform_seq,
        "mrna_seq": mrna_seq,
        "acid": acids,
        "exon_flag": format_flag(exons_flag),
        "cds_flag": format_flag(cds_flag),
        "junc_flag_for_acid": format_flag(junc_flag),
    }
    for key, value in full_info.items():
        df_9mer[key] = value
    return df_9mer


def preprocess(cds_gff, junction_file, dna_seq, protein_seq):
    position_info = get_position_info(cds_gff)

    junctions = get_novel_junction(junction_file)
    novel_isoform_ids = junctions["isoform"].unique()
    print(novel_isoform_ids.shape)

    isoform = SeqIO.parse(dna_seq, "fasta")
    isoform = pd.DataFrame(
        data=[(record.id, str(record.seq)) for record in isoform],
        columns=["transcript_id", "seq"],
    )
    isoform.set_index("transcript_id", inplace=True)
    assert set(novel_isoform_ids).issubset(set(isoform.index))

    isoform = isoform.loc[novel_isoform_ids]
    isoform["seq"] = isoform["seq"].str.upper()

    peptides = SeqIO.parse(protein_seq, "fasta")
    peptides = pd.DataFrame(
        data=[(record.id, str(record.seq)) for record in peptides],
        columns=["transcript_id", "peptide"],
    )
    peptides["peptide"] = peptides["peptide"].str.upper()

    print(
        f"{protein_seq.split('/')[-1]} contains {len(peptides)} rows, but only {len(peptides['transcript_id'].unique())} unique transcript ids"
    )
    print(
        f"these lost in isoform: {pd.DataFrame(set(peptides['transcript_id']) - set(isoform.index))}"
    )
    print(
        f"these duplicated in pepeptides: {peptides[peptides.duplicated(subset='transcript_id', keep=False)]}"
    )
    peptides.set_index("transcript_id", inplace=True)

    all_transcript_ids = list(set(isoform.index) & set(peptides.index))
    if len(all_transcript_ids) == 0:
        raise ValueError("no novel isoform has peptides")
    print(f"In total {len(isoform)} novel isoforms, and total {len(peptides)} peptides")
    print(f"{len(all_transcript_ids)} novel isoforms have paired peptides")
    return all_transcript_ids, position_info, isoform, peptides, junctions


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate novel junction derived 9-mer neoantigens"
    )
    parser.add_argument(
        "--cds_gff",
        default="example/example.cds.gtf",
        help="Path to CDS GFF file",
    )
    parser.add_argument(
        "--junction_file",
        default="example/example_junctions.txt",
        help="Path to NPC_junctions.txt",
    )
    parser.add_argument(
        "--dna_seq",
        default="example/example.fasta",
        help="Path to corrected DNA FASTA",
    )
    parser.add_argument(
        "--protein_seq",
        default="example/example.faa",
        help="Path to protein FASTA",
    )
    parser.add_argument(
        "--output_file",
        default="./output/positive_novel_junction_peptides_example.csv",
        help="Path to output CSV",
    )
    parser.add_argument(
        "--parallel",
        action="store_true",
        default=False,
        help="Enable parallel processing (default: False)",
    )
    args = parser.parse_args()

    all_transcript_ids, position_info, isoform, peptides, junctions = preprocess(
        args.cds_gff, args.junction_file, args.dna_seq, args.protein_seq
    )

    results = []
    if not args.parallel:
        for transcript_id in all_transcript_ids[:]:
            exons_info = position_info[
                (position_info["transcript_id"] == transcript_id)
                & (position_info["type"] == "exon")
            ].copy()
            cds_info = position_info[
                (position_info["transcript_id"] == transcript_id)
                & (position_info["type"] == "CDS")
            ].copy()
            isoform_seq = isoform.loc[transcript_id, "seq"]
            acids_seq = peptides.loc[transcript_id]["peptide"]
            if not isinstance(acids_seq, str):
                acids_seqs = acids_seq.tolist()
            else:
                acids_seqs = [acids_seq]

            success_times = 0
            for acids_seq in acids_seqs:
                junction_id = junctions.loc[
                    junctions["isoform"] == transcript_id, "junction_number"
                ]
                junction_id = junction_id.apply(lambda x: int(x.split("_")[1])).values
                df_9mer = run_one_transcript(
                    transcript_id,
                    junction_id,
                    isoform_seq,
                    acids_seq,
                    exons_info,
                    cds_info,
                )
                if df_9mer is not None:
                    results.append(df_9mer)
                    success_times += 1
            if success_times == 0:
                raise ValueError(f"transcript_id {transcript_id} has no valid peptide")
    else:
        process_pool = pathos.pools.ProcessPool(150)
        for transcript_id in all_transcript_ids[:]:
            exons_info = position_info[
                (position_info["transcript_id"] == transcript_id)
                & (position_info["type"] == "exon")
            ].copy()
            cds_info = position_info[
                (position_info["transcript_id"] == transcript_id)
                & (position_info["type"] == "CDS")
            ].copy()
            isoform_seq = isoform.loc[transcript_id, "seq"]
            true_acids = peptides.loc[transcript_id]["peptide"]
            junction_id = junctions.loc[
                junctions["isoform"] == transcript_id, "junction_number"
            ]
            junction_id = junction_id.apply(lambda x: int(x.split("_")[1])).values
            async_result = process_pool.apipe(
                run_one_transcript,
                transcript_id,
                junction_id,
                isoform_seq,
                true_acids,
                exons_info,
                cds_info,
            )
            results.append(async_result)
        for index in range(len(results)):
            results[index] = results[index].get()
        process_pool.close()
        process_pool.join()

    results = pd.concat(results, axis=0, ignore_index=True)
    Path(args.output_file).parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.output_file, index=False, sep="\t")
