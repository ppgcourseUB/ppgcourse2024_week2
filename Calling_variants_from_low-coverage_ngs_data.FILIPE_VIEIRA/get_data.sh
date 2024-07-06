#!/bin/bash

export REG="11:20000000-23000000"

samtools faidx http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz 11 > hs37d5.fa.gz
samtools faidx hs37d5.fa.gz

#parallel --header : "samtools view -o {sample}.{pop}.bam http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/{1}/alignment/{sample}.chrom11.ILLUMINA.bwa.{pop}.low_coverage.{date}.bam $REG" :::: samples.annot
tail -n +2 samples.annot | xargs -l -P 20 bash -c 'samtools view -o $0.$1.bam http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/$0/alignment/$0.chrom11.ILLUMINA.bwa.$1.low_coverage.$4.bam $REG'
rm *.bai
#parallel samtools index ::: *.bam
ls *.bam | xargs -P 20 -n 1 samtools index

ls *.bam | sort -t "." -k 2,2 > samples.bam_list
