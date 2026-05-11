#!/bin/bash


input_dir=/path/to/raw_data

sample_list=/path/to/sample_list.txt

output_dir=/path/to/output/optitype

hla_fasta=/path/to/hla_reference_rna.fasta

cat ${sample_list} | while read id
do
    mkdir -p ${output_dir}/${id}/

    if [ ! -s ${output_dir}/${id}/${id}_1.hla.bam ]; then
        echo "Running razers3 on ${id}_1.fq.gz"
        razers3 -i 95 -m 1 -dr 0 \
            -o ${output_dir}/${id}/${id}_1.hla.bam \
            ${hla_fasta} \
            ${input_dir}/${id}_1.fq.gz \
            1>> ${output_dir}/${id}/${id}_optitype.log 2>&1
    fi

    if [ ! -s ${output_dir}/${id}/${id}_2.hla.bam ]; then
        echo "Running razers3 on ${id}_2.fq.gz"
        razers3 -i 95 -m 1 -dr 0 \
            -o ${output_dir}/${id}/${id}_2.hla.bam \
            ${hla_fasta} \
            ${input_dir}/${id}_2.fq.gz \
            1>> ${output_dir}/${id}/${id}_optitype.log 2>&1
    fi

samtools bam2fq -@ 16 ${output_dir}/${id}/${id}_1.hla.bam > ${output_dir}/${id}/${id}_1.hla.fq
samtools bam2fq -@ 16 ${output_dir}/${id}/${id}_2.hla.bam > ${output_dir}/${id}/${id}_2.hla.fq

rm ${output_dir}/${id}/${id}_1.hla.bam
rm ${output_dir}/${id}/${id}_2.hla.bam

mkdir -p ${output_dir}/${id}/optitype/
echo "Running OptiType on ${id}"
OptiTypePipeline.py \
    -i ${output_dir}/${id}/${id}_1.hla.fq ${output_dir}/${id}/${id}_2.hla.fq \
    --rna -v \
    -o ${output_dir}/${id}/optitype/

done


