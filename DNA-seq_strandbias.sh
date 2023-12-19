#!usr/bin/env bash
conda activate conda-env
#bwa index $ref

sample='Y1'
raw_reads='/storage/coda1/p-fstorici3/0/shared/raw_reads/DNA-seq/FS41/' #keep your fq files here
trimmed_reads='trimmed'
aligned='aligned'
vcf='vcf'
DB='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/omics-data/DNA-seq/FSlab/DBs/'
mkdir $trimmed_reads $aligned $aligned/log $aligned/fixmate $aligned/temp


run() {
####################Trimming####################
    #trim_galore -q 15 --length 50 --fastqc --paired ${raw_reads}/${sample}_R1.fq ${raw_reads}/${sample}_R2.fq -o $trimmed_reads 
    #bwa mem -R '@RG\tID:'${sample}'\tSM:bar\tLB:1X' $ref ${trimmed_reads}/${sample}_R1_val_1.fq ${trimmed_reads}/${sample}_R2_val_2.fq > ${aligned}/${sample}.sam
    bowtie2 --threads 4 -x $(dirname $ref)/$(basename $ref .fa) -1 ${trimmed_reads}/${sample}_R1_val_1.fq -2 ${trimmed_reads}/${sample}_R2_val_2.fq -S ${aligned}/${sample}.sam 2> ${aligned}/log/${sample}.log
    #samtools fixmate -O bam ${aligned}/${sample}.sam ${aligned}/fixmate/fix_${sample}.bam
    samtools sort -O bam -T $aligned/temp/ ${aligned}/fixmate/fix_${sample}.bam | samtools addreplacerg -r '@RG\tID:samplename\tSM:samplename' - -o ${aligned}/${sample}.bam
    #samtools view -b ${aligned}/fixmate/fix_${sample}.bam | samtools sort -O bam -@ 4 - -o ${aligned}/${sample}.bam
    #samtools addreplacerg -r '@RG\tID:samplename\tSM:samplename' ${aligned}/${sample}.bam -o ${aligned}/final_${sample}.bam
    samtools index ${aligned}/${sample}.bam
    #gatk CreateSequenceDictionary -R $ref
    gatk --java-options "-Xmx8G" HaplotypeCaller -R $ref -I ${aligned}/${sample}.bam -O ${vcf}/haplo_${sample}.vcf.gz #Germline variants
    gatk Mutect2 -R $ref -I ${aligned}/${sample}.bam --germline-resource ${DB}/af-only-gnomad.hg38.vcf.gz --panel-of-normals ${DB}/1000g_pon.hg38.vcf.gz -O ${vcf}/mutect_${sample}.vcf.gz #Somatic variants for humans
    #gatk Mutect2 -R $ref -I ${aligned}/${sample}.bam -O ${vcf}/mutect_${sample}.vcf.gz #All mutations ofr yeast tumor mode only

    #gatk Mutect2 -R $ref -I ${aligned}/${sample}.bam --germline-resource af-only-gnomad.hg38.vcf.gz --panel-of-normals 1000g_pon.hg38.vcf.gz -O ${vcf}/${sample}.vcf.gz
}

ref=$hgref
for sample in CD4T-1 CD4T-2 CD4T-3 CD4T-4 DLTB-8 DLTB-P TLTB-8 TLTB-P H9-1 H9-2 HEK293T-RNASEH2A-KO-T17 HEK293T-RNASEH2A-KO-T8 HEK293T-WT 
do 
run &
done

ref=$saccerref
for sample in E134-top1-KO E134-top1-rnh201-KO E134-WT BY4742-top1-KO BY4742-top1-rnh201-KO BY4742-WT
do
run &
done

ref=$saccerref
for sample in Y1 Y2 Y3 Y4 Y5 Y6
do
run &
done

/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/bin/gatk-4.3.0.0

ref=$hgref
for sample in H1 H2 H3 H4 H5
do
bowtie2 --threads 4 -x $(dirname $ref)/$(basename $ref .fa) -1 ${trimmed_reads}/${sample}_R1_val_1.fq -2 ${trimmed_reads}/${sample}_R2_val_2.fq -S ${aligned}/${sample}.sam 2> ${aligned}/log/${sample}.log &
done

####################Strand bias info####################
strand_bias() {
    forward=$(samtools view $file chrM | gawk '(and(16,$2))' | wc -l)
    reverse=$(samtools view $file chrM | gawk '(! and(16,$2))' | wc -l)
    f=$(basename $file .bam)
    echo $f'\t'$forward'\t'$reverse
}

strand_bias_all() {
    forward=$(samtools view $file | gawk '(and(16,$2))' | wc -l)
    reverse=$(samtools view $file | gawk '(! and(16,$2))' | wc -l)
    f=$(basename $file .bam)
    echo $f $forward $reverse
}

genomiccoverage() {
    samtools depth -r chrM:1-16569 $file | awk -v OFS="\t" -F"\t" '{sum+=$3} END { print $file ,sum/NR}'
}

for file in $(ls aligned/*.bam)
do
strand_bias 
done

cd aligned
for file in $(ls HEK*.bam)
do
samtools view -b $file "chr19:12807584-12814640" > rnh2_$(basename $file)
samtools fasta rnh2_$(basename $file) > rnh2_$(basename $file .bam).fasta
done
