#!/bin/bash

input_gtf=/path/to/merged.combined.adjusted.gtf

ref_gtf=/path/to/gencode.annotation.gtf

genome_fasta=/path/to/genome.fa

cage_peak=/path/to/refTSS_v3.3_human_coordinate.hg38.bed

polyA_motif=/path/to/mouse_and_human.polyA_motif.txt

polyA_peak=/path/to/atlas.clusters.2.0.GRCh38.96.bed

coverage=/path/to/SJ.out.tab

prefix=project

qc_dir=/path/to/output/SQANTI3_QC

sqanti3_qc.py -t 30 ${input_gtf} ${ref_gtf} ${genome_fasta} \
    --CAGE_peak ${cage_peak} \
    --polyA_motif_list ${polyA_motif} \
    --polyA_peak ${polyA_peak} \
    -c ${coverage} \
    -o ${prefix} \
    -d ${qc_dir} \
    --cpus 4 \
    --report both
