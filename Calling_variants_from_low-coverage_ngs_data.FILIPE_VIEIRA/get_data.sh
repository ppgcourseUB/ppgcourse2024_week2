#!/bin/bash

# Download reference (chr11)
samtools faidx http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz 11 > hs37d5.fa.gz
samtools faidx hs37d5.fa.gz

# Download example BAM files (chr11, 1000G)
#parallel --header : "wget -O {sample}.{pop}.bam http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/{1}/alignment/{sample}.chrom11.ILLUMINA.bwa.{pop}.low_coverage.{date}.bam" :::: samples.annot
#parallel samtools index ::: *.bam

tail -n +2 samples.annot | xargs -l -P 20 bash -c 'wget -O $0.$1.bam http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/$0/alignment/$0.chrom11.ILLUMINA.bwa.$1.low_coverage.$4.bam'
ls *.bam | xargs -P 20 -n 1 samtools index

# Create BAM list
ls *.bam | sort -t "." -k 2,2 > samples.bam_list
