#!/usr/bin/env python

"""
Please change the path for each program
12/02/2019
Veronica Mixao
CRG (Barcelona)

from programs_config import program_path
"""

#Please set the correct path for each program
path = {}

path["haplotypo"] = "~/src/haplotypo/bin/"
path["java"] = "java" 
path["bwa"] = "bwa" 
path["samtools"] = "samtools" 
path["picard"] = "~/src/haplotypo/dependencies/gatk-4.0.12.0/gatk"
path["gatk"] = "~/src/haplotypo/dependencies/gatk-4.0.12.0/gatk-package-4.0.12.0-local.jar" 
path["freebayes"] = "freebayes "
path["bcftools"] = "bcftools" 
path["vcffilter"] = "vcffilter"

def program_path(program):
	return path[program]

#haplotypo was tested with the following program versions
#java v1.8.0_171
#bwa v0.7.15
#samtools v1.6
#picard v4.0.2.1
#gatk v4.0.2.1
#freebayes v1.1.0-50-g61527c5
#bcftools v1.9
