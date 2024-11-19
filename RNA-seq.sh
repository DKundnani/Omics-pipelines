#!usr/bin/bash

################### Input from the user
raw_reads='/storage/coda1/p-fstorici3/0/shared/raw_reads/RNA-seq/24062-07/fqs' #Directory where fq data is present
saf=$(dirname $raw_reads)/featurecounts/exons.saf #tsv file with annotations: Cols as: GeneID,chr,start,stop,strand
ref=~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa #fasta file of reference genome
samplelist='HEK293T-WT-1 HEK293T-WT-2 HEK293T-KO-T3-8-1 HEK293T-KO-T3-8-2 HEK293T-KO-T3-17-1 HEK293T-KO-T3-17-2'

######### Ouput folders created
trimmed=$(dirname $raw_reads)/trimmed; mkdir $trimmed
bam=$(dirname $raw_reads)/bams; mkdir $bam; mkdir $bam/log; mkdir $bam/fixmate; mkdir $bam/temp/
featurecounts=$(dirname $raw_reads)/featurecounts ; mkdir $featurecounts


rnaseq () {
    #1. Trimming using TrimGalore
    trim_galore -q 15 --length 50 --fastqc --paired ${raw_reads}/${sample}_R1.fq ${raw_reads}/${sample}_R2.fq -o $trimmed; echo $sample Trimming Done!
    #2. Alignment using Bowtie2
    bowtie2 --threads 16 -x $(dirname $ref)/$(basename $ref .fa) -1 ${trimmed}/${sample}_R1_val_1.fq -2 ${trimmed}/${sample}_R2_val_2.fq -S ${bam}/${sample}.sam 2> ${bam}/log/${sample}.log; echo $sample Alignment Done!
    #grep overall $bam/log/${sample}.log
    #3. Assign pairs for paired end reads using samtools fixmate
    samtools fixmate -O bam ${bam}/${sample}.sam ${bam}/fixmate/fix_${sample}.bam; echo $sample Fixmate Done!
    #4. Adding sample name, Sorting #samtools addreplace gives sam as an output so keep it first and then sort.
    samtools addreplacerg -r $(echo @RG\\tID:${sample}\\tSM:${sample}) ${bam}/fixmate/fix_${sample}.bam | samtools sort -O bam -T $bam/temp/ - -o ${bam}/${sample}.bam
    samtools sort -O bam -T $bam/temp/ ${bam}/fixmate/fix_${sample}.bam -o ${bam}/${sample}.bam; echo $sample Resorting Done!
    #5. Indexing
    samtools index ${bam}/${sample}.bam; echo $sample Indexing Done!
}



for sample in $samplelist; do rnaseq & done
wait
grep overall $bam/log/*
rm -r ${bam}/fixmate ${bam}/temp
#5. Counting RNA-seq reads in gene locations in sacCer3 genome using featureCounts
featureCounts -p -O -T 8 -F SAF -a $saf -o $featurecounts/$(basename $saf .saf).counts $bam/*bam


