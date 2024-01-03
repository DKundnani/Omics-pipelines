#!usr/bin/env python3

#Usage: Create a vcf folder to plaice all vcf files and run from outside that folder. 

from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
import sigProfilerPlotting as sigPlt

genInstall.install('yeast', rsync=False, bash=True)

matrices = matGen.SigProfilerMatrixGeneratorFunc("vcf", "yeast", "./vcf/",plot=True, exome=False, bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)

sigPlt.plotID("vcf/output/ID/vcf.ID83.all", "vcf/output/", "vcf", "83", percentage=True)

