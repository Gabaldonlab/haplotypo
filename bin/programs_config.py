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

path["haplotypo"] = "~/users/tg/vdepinho/Candida_albicans/technical_analysis/Xpipeline/"
path["java"] = "/usr/bin/java" 
path["bwa"] = "~/users/tg/vdepinho/Programs/bwa/bwa-0.7.15/bwa" 
path["samtools"] = "~/users/tg/vdepinho/Programs/samtools/version1.9/samtools-1.9/samtools" 
path["picard"] = "~/users/tg/vdepinho/Programs/GATK/GATK_v4.01/gatk-4.0.2.1/gatk "
path["gatk"] = "~/users/tg/vdepinho/Programs/GATK/GATK_v4.01/gatk-4.0.2.1/gatk-package-4.0.2.1-local.jar" 
path["freebayes"] = "~/users/tg/lcarrete/src/freebayes/bin/freebayes "
path["bcftools"] = "~/users/tg/vdepinho/Programs/samtools/version1.9/bcftools-1.9/bcftools" 
path["vcffilter"] = "~/users/tg/cpegueroles/soft/vcflib/bin/vcffilter"

	
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
