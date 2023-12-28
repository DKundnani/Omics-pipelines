#source ~/.bash_profile
#Initializing variables
in=/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/omics-data/WGBS_CD4T #folder containing fastq reads
trimmed_reads=/storage/home/hcoda1/5/dkundnani3/scratch/Bismark/trimmed #folder for trimmed reads
samplelist='WT K46W G42S rnh201-null'
### trimming PE reads
for sample in $samplelist; do
R1=${in}/${sample}/${sample}_R1.fastq.gz
R2=${in}/${sample}/${sample}_R2.fastq.gz
#conda activate conda-env #environment with trim-galore
trim_galore -q 15 --length 50 --fastqc --paired $R1 $R2 -o $trimmed_reads &
done
wait
