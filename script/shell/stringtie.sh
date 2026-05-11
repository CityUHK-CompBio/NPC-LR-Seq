#!/bin/bash

sample_list=/path/to/sample_list.txt

reference_gtf=/path/to/reference/gencode38.annotation.gtf

bam_dir=/path/to/star_alignment/uniquely_map

stringtie_out=/path/to/output/StringTie
mkdir -p ${stringtie_out}

cat ${sample_list} | while read sample
do
    input_file=${bam_dir}/${sample}Aligned.sortedByCoord.out.bam

    if [ ! -s ${stringtie_out}/${sample}.gtf ]; then
        stringtie -G ${reference_gtf} -p 16 \
            -o ${stringtie_out}/${sample}.gtf \
            ${input_file}
    fi
done



stringtie_out=/path/to/output/StringTie

reference_gtf=/path/to/reference/gencode38.annotation.gtf

merged_gtf=/path/to/output/StringTie/merged.gtf

ls "${stringtie_out}"/*.gtf | grep -v '/merged\.gtf$' > gtf_list.txt

stringtie --merge -p 10 \
    -G ${reference_gtf} \
    -o ${merged_gtf} \
    -l merge gtf_list.txt 
