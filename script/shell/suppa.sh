#!/bin/bash

gtf=/path/to/output/gtf/merged.gtf

out_prefix=/path/to/output/suppa/project.suppa

suppa.py generateEvents \
    -i ${gtf} \
    -o ${out_prefix} \
    -f ioe \
    -e SE AF MX RI A5 A3 AL



awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' /path/to/output/suppa/*.ioe \
> /path/to/output/suppa/project.all.events.ioe


ioe_merge_file=/path/to/output/suppa/project.all.events.ioe

expression_matrix=/path/to/expression/isoform-count-matrix.txt

out_prefix=/path/to/output/suppa/project_events

log_file=/path/to/output/suppa/psiPerEvent_log.txt

suppa.py psiPerEvent \
    --ioe ${ioe_merge_file} \
    --expression-file ${expression_matrix} \
    -o ${out_prefix} \
    1> ${log_file} 2>&1 



ioe_merge_file=/path/to/output/suppa/project.all.events.ioe

psi_cond1=/path/to/condition1.psi
psi_cond2=/path/to/condition2.psi

tpm_cond1=/path/to/condition1.tpm
tpm_cond2=/path/to/condition2.tpm

out_prefix=/path/to/output/suppa/project_diffSplice

suppa.py diffSplice \
    -m empirical -gc \
    --input ${ioe_merge_file} \
    --save_tpm_events \
    -p ${psi_cond1} ${psi_cond2} \
    --tpm ${tpm_cond1} ${tpm_cond2} \
    -o ${out_prefix}
