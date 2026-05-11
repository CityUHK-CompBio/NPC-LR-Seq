#!/bin/bash

run_star_uniquely(){
cat $1 | while read id                                                              
do
mkdir -p $3
read1=$2/${id}_1.fastq.gz                                 
read2=$2/${id}_2.fastq.gz                                                                         
STARgenomeDir=$4                                                                                       
nThreadsSTAR=16                                                                                                       

if grep -q ${id} $2/unpaired.txt;then
STARparCommon="--genomeDir $STARgenomeDir  --readFilesIn $read1  --outSAMunmapped Within \
--outFilterType BySJout --outSAMattributes Standard   --outFilterMultimapNmax 1 --outSAMattrIHstart 1  \ 
--outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04  --outSAMattrRGline ID:${id} SM:${id} PL:ILLUMINA \
--alignIntronMin 20 --alignIntronMax 1000000   --alignMatesGapMax 1000000 --alignSJoverhangMin 8 \  
--alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesCommand zcat --outSAMstrandField intronMotif "

else
STARparCommon="--genomeDir $STARgenomeDir  --readFilesIn $read1 $read2   --outSAMunmapped Within \
--outFilterType BySJout --outSAMattributes  Standard   \ 
--outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04  --outSAMattrRGline ID:${id} SM:${id} PL:ILLUMINA \
--alignIntronMin 20 --alignIntronMax 1000000   --alignMatesGapMax 1000000 --alignSJoverhangMin 8 \  
--alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesCommand zcat --outSAMstrandField intronMotif --outFilterMultimapNmax 1 --outSAMattrIHstart 1 "                           

fi
echo $STARparCommon
if [ ! -s $3/${id}Log.final.out ];then                                                                   
STARparRun=" --runThreadN $nThreadsSTAR --genomeLoad LoadAndKeep"                      
STARparBAM="--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFileNamePrefix $3/${id}"                    
STARparsMeta="--outSAMheaderCommentFile commentsENCODEli:ong.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate --limitBAMsortRAM 411432265264" 
                                                                                                                      
echo ${id}                                                                                                   
STAR $STARparCommon $STARparRun $STARparBAM $STARparStrand $STARparsMeta >> $3/star.log 2>&1                                              
fi

done  
}



