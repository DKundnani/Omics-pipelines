#Getting consesnsus

ref='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/reference/sacCer2/sacCer2.fa'
vcfpath='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/AGS/DNAseq/vcf/vcf'
consensus_output='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/AGS/DNAseq/consensus'
coordinates=$consensus_output/coordinates

for sample in Y1 Y2 Y5 Y6; do
cat $ref | bcftools consensus ${vcfpath}/haplo_${sample}.vcf.gz > ${consensus_output}/${sample}.fa &
done

## coordiantes has co-ordinates of specific region of genome
## chrII    4470    5511

for sample in Y1 Y2 Y5 Y6; do
bedtools getfasta -fi ${consensus_output}/${sample}.fa -bed ${coordinates} -fo ${consensus_output}/${sample}_regions.fa &
done
