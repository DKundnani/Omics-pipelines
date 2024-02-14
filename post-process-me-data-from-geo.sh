#!usr/bin/env bash

conda activate conda-env

cut -f1-3,6 HEK293T_WGBS_hg38.bed | grep -Ev 'random|alt|Un|M|X|Y' | sort -k1,1 -k2,2n | uniq >  HEK293T_WGBS_hg38.bg
#lift over if necessary

## Get unique locations of methylation
conda activate r_env
## In R
conda activate r_env 
R #Open R terminal
df=read.table('HEK293T_WGBS_hg38.bg', sep='\t', header=FALSE)
library(dplyr)
#mydf_2 <- df %>% group_by(color) %>% summarize(value=mean(value))
finaldf=aggregate(V4 ~ V1 + V2 +V3, data = df, FUN = mean, na.rm = TRUE)
write.table(finaldf, 'HEK293T_WGBS_hg38_unique.bg' , sep = "\t",quote=F,row.names=F, col.names=F)

## final bw file
conda activate conda-env
bedGraphToBigWig HEK293T_WGBS_hg38_unique.bg $hggenome HEK293T_WGBS_hg38_unique.bw
