#Initializing variables
conda activate conda-env #activate conda environment with bismark tool 
ref='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa' #reference fasta file
aligned='/storage/home/hcoda1/5/dkundnani3/scratch/Bismark/bams' #folder for alignment files
trimmed_reads=''
###################################################
# Start pipeline
###################################################
tR1=$trimmed_reads/$(basename $R1 .fq.gz)_val_
tR2=$trimmed_reads/
pigz -dfkp 6 $FQZ1 & pigz -dfkp 6 $FQZ2
FQ1=`echo $FQZ1 | sed 's/.gz$//'`
FQ2=`echo $FQZ2 | sed 's/.gz$//'`
wait

###################################################
# Run Bismark with bowtie2 alignment
###################################################
bismark --non_directional --bowtie2 -gzip \
-p 4 $REF -1 $FQ1 -2 $FQ2
rm $FQ1 $FQ2
rm ${FQ1}_C_to_T.fastq ${FQ2}_G_to_A.fastq
pigz ${FQ1}_bismark_bt2_pe.sam
}

export -f wgbsaln
parallel -u -j2 wgbsaln ::: *pair1.fastq.gz

#Jobs can also be run over several servers with the following setup (further reading here).
#parallel -u -j2 -S server1,server2,server3 wgbsaln ::: *pair1.fastq.gz

#The alignment will take quite a bit longer than a standard gDNA alignment. When it's finished, the compressed sam files can be concatenated (without uncompressing):

cat Sample1_*.sam.gz > Sample1.sam.gz

#These compressed sam files can be used directly in methylation calling, but before we start, take a moment to look at the M-bias graphs, as suggested in the manual to determine whether the first and last few bases of read 1 and 2 show any strong biases. If so, you'll need to ignore those bases in the methylation calling. You might also want to move the concatenated sam.gz alignments to a new directory to perform methylation calling with the following script.

#################################################
# Extract methylation scores
#################################################
methextr(){
REF=/pathto/refgenome/
SAMGZ=$1
SAM=`echo $SAMGZ | sed 's/.gz//'`
bismark_methylation_extractor -p --bedgraph --counts \
 --buffer_size 20G --cytosine_report --genome_folder \
$REF --multicore 10 --gzip --ignore 3 --ignore_r2 3 $SAMGZ
}
export -f methextr
parallel -u -j3 methextr ::: Sample1.sam.gz Sample2.sam.gz SampleN.sam.gz

#The script is based upon an example in the manual that generates a genome wide CpG report in addition to the Bedgraph output that look like this.
'''
==> Sample1.sam.bismark.cov <==
20 60179 60179 34.4827586206897 10 19
20 60180 60180 46.6666666666667 7 8
20 60358 60358 0 0 1
20 60426 60426 92.3076923076923 12 1
20 60427 60427 92 23 2
20 60432 60432 84.6153846153846 11 2
20 60433 60433 72 18 7
20 60551 60551 38.8888888888889 7 11
20 60552 60552 46.6666666666667 14 16
20 60578 60578 47.3684210526316 9 10

==> Sample1.sam.CpG_report.txt <==
20 60179 + 10 19 CG CGT
20 60180 - 7 8 CG CGT
20 60426 + 12 1 CG CGA
20 60427 - 23 2 CG CGG
20 60432 + 11 2 CG CGA
20 60433 - 18 7 CG CGG
20 60551 + 7 11 CG CGA
20 60552 - 14 16 CG CGT
20 60578 + 9 10 CG CGT
20 60579 - 25 15 CG CGC
'''
   

#Reference: http://genomespot.blogspot.com/2015/07/genome-methylation-analysis-with-bismark.html; https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf; https://github.com/FelixKrueger/Bismark 

