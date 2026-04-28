#!/bin/bash

index_build() {
    mkdir -p "$3"
    rsem-prepare-reference --gtf "$1" "$2" "$3/$4" -p 8
}

quantify() {
    local sample_list="$1"
    local align_dir="$2"
    local rsem_index="$3"
    local out_dir="$4"
    local unpaired_file

    mkdir -p "$out_dir"
    unpaired_file="$(dirname "$sample_list")/unpaired.txt"

    while IFS= read -r sample; do
        [ -z "$sample" ] && continue

        local result_file="$out_dir/${sample}.isoforms.results"
        local bam_file="$align_dir/${sample}Aligned.toTranscriptome.out.bam"

        [ -s "$result_file" ] && continue
        [ -s "$bam_file" ] || {
            echo "Warning: missing BAM for ${sample}, skip"
            continue
        }

        echo "Quantifying $sample"

        if [ -s "$unpaired_file" ] && grep -Fxq "$sample" "$unpaired_file"; then
            echo "${sample} single-end"
            rsem-calculate-expression --no-bam-output \
                --strandedness none \
                --forward-prob 0.5 \
                --alignments \
                -p 8 \
                "$bam_file" \
                "$rsem_index" \
                "$out_dir/$sample"
        else
            echo "${sample} paired-end"
            rsem-calculate-expression --paired-end \
                --no-bam-output \
                --strandedness none \
                --forward-prob 0.5 \
                --alignments \
                -p 8 \
                "$bam_file" \
                "$rsem_index" \
                "$out_dir/$sample"
        fi
    done < "$sample_list"
}


