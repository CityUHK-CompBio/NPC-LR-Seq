#!/bin/bash

lr_dirs=/path/to/LR_seq_data/

sr_dirs=/path/to/SR_seq_data/

lr_samples=/path/to/lr_sample.txt

sr_samples=/path/to/sr_sample.txt

FA=/path/to/Iso_Seq_Primer.fa

Genome=/path/to/genome.fa

iso_preprocess(){
    source /path/to/miniconda3/bin/activate isoseq3_env

    outdir=$1
    mkdir -p ${outdir}

    sample_list=/path/to/lr_samples.txt

    cat ${sample_list} | while read sample
    do
        if [ -s ${outdir}/${sample}.clustered.hq.fasta.gz ]; then
            echo "${sample} has already been processed"
        else
            echo "Processing ${sample}"

            if [ ! -s ${outdir}/${sample}.ccs.bam ]; then 
                echo "Generating CCS reads for ${sample}"
                ccs ${lr_dirs}/${sample}.subreads.bam \
                    ${outdir}/${sample}.ccs.bam \
                    --min-rq 0.9 -j 64
            fi

            if [ -s ${outdir}/${sample}.ccs.bam ] && [ ! -s ${outdir}/${sample}.flnc.bam ]; then
                echo "Removing primers and identifying barcodes for ${sample}"
                lima ${outdir}/${sample}.ccs.bam \
                     ${FA} \
                     ${outdir}/${sample}.fl.bam \
                     --isoseq --peek-guess -j 64

                echo "Refining ${sample}"
                isoseq3 refine ${outdir}/${sample}.*NEB_Clontech_3p.bam \
                               ${FA} \
                               ${outdir}/${sample}.flnc.bam -j 64
            fi

            if [ -s ${outdir}/${sample}.flnc.bam ] && [ ! -s ${outdir}/${sample}.clustered.bam ]; then
                echo "Clustering ${sample}"
                isoseq3 cluster ${outdir}/${sample}.flnc.bam \
                                ${outdir}/${sample}.clustered.bam \
                                --verbose --use-qvs --poa-cov 100 -j 64
            fi
        fi
    done
}

lr_dirs=/path/to/LR_seq_data
FA=/path/to/Iso_Seq_Primer.fa

outdir=/path/to/output/iso_preprocess
iso_preprocess ${outdir}


alignment(){
    indir=$1
    outdir=$2
    mkdir -p ${outdir}

    source /path/to/miniconda3/bin/activate isoseq3_env

    sample_list=/path/to/lr_samples.txt

    cat ${sample_list} | while read sample
    do
        if [ -s ${indir}/${sample}.clustered.hq.bam ] && [ ! -s ${outdir}/${sample}.clustered.hq.sorted.bam ]; then
            echo "Aligning ${sample} clustered HQ reads to reference genome"

            pbmm2 align --preset ISOSEQ --sort \
                ${Genome} \
                ${indir}/${sample}.clustered.hq.bam \
                ${outdir}/${sample}.clustered.hq.sorted.bam
        fi
    done
}

indir=/path/to/output/iso_preprocess
outdir=/path/to/output/pbmm2_alignment

Genome=/path/to/genome.fa

alignment ${indir} ${outdir}

collapse(){
    indir=$1
    outdir=$2
    mkdir -p ${outdir}

    source /path/to/miniconda3/bin/activate isoseq3_env

    sample_list=/path/to/lr_samples.txt

    cat ${sample_list} | while read sample
    do
        if [ -s ${indir}/${sample}.clustered.hq.sorted.bam ] && [ ! -s ${outdir}/${sample}.collapse.gff ]; then
            echo "Collapsing isoforms for ${sample}"
            isoseq3 collapse \
                ${indir}/${sample}.clustered.hq.sorted.bam \
                ${outdir}/${sample}.collapse.gff \
                --min-aln-coverage 0.99 \
                --min-aln-identity 0.95
        fi
    done
}

indir=/path/to/output/pbmm2_alignment
outdir=/path/to/output/collapse

collapse ${indir} ${outdir}



sample_list=/path/to/lr_samples.txt
collapse_dir=/path/to/output/collapse

out_file=${collapse_dir}/num_full_transcript.txt

echo -e "sample\tnumber" > ${out_file}

cat ${sample_list} | while read id
do
    if [ -s ${collapse_dir}/${id}.collapse.gff ]; then
        num=$(grep -w "transcript" ${collapse_dir}/${id}.collapse.gff | wc -l)
        echo -e "${id}\t${num}" >> ${out_file}
    else
        echo "Warning: ${collapse_dir}/${id}.collapse.gff not found, skipping"
    fi
done

gffcompare(){
    source /path/to/miniconda3/bin/activate gffcompare_env

    indir=$1
    outdir=$2
    mkdir -p ${outdir}

    ls ${indir}/*.collapse.gff > ${outdir}/input_list

    reference_gtf=/path/to/gencode.annotation.gtf

    command gffcompare \
        -r ${reference_gtf} \
        -M -S \
        -o ${outdir}/merged \
        -i ${outdir}/input_list
}

indir=/path/to/output/collapse
outdir=/path/to/output/gffcompare

gffcompare ${indir} ${outdir}


sqanti3_qc(){
    source /path/to/miniconda3/bin/activate sqanti3_env

    merged_gtf=$1       # e.g. gffcompare merged.combined.gtf
    reference_gtf=$2    # reference annotation (e.g. GENCODE)
    reference_fasta=$3  # reference genome fasta
    out_prefix=$4       # prefix for output files
    out_dir=$5          # output directory

    mkdir -p ${out_dir}

    sqanti3_qc.py \
        -t 30 \
        ${merged_gtf} \
        ${reference_gtf} \
        ${reference_fasta} \
        -o ${out_prefix} \
        -d ${out_dir} \
        --cpus 4 \
        --report both
}

merged_gtf=/path/to/output/gffcompare/merged.combined.gtf
reference_gtf=/path/to/reference/gencode.annotation.gtf
reference_fasta=/path/to/reference/genome.fa
out_prefix=project_name
out_dir=/path/to/output/sqanti3

sqanti3_qc ${merged_gtf} ${reference_gtf} ${reference_fasta} ${out_prefix} ${out_dir}
