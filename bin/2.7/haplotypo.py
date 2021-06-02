#!/usr/bin/env python

"""
Pipeline: Read mapping and variant calling on phased genomes with correction for variants and respective haplotype reconstruction
12/02/2018
Veronica Mixao 
CRG (Barcelona)
"""

import argparse
import os
import sys
import configparser
from programs_config import program_path

#Setting the correct path for each program
haplotypo = program_path("haplotypo")
if program_path("python") == None:
	sys.exit("ERROR!! Please revise the Python path parameters!")
python_path = program_path("python")

##############
#### MAIN ####
##############

parser = argparse.ArgumentParser(description="Read mapping and variant calling on phased genomes with correction for variants and haplotype reconstruction.")
parser.add_argument("-o", "--outdir", dest="outDir", action="store", help="Directory where the data will be stored")
parser.add_argument("-thr", "--threads", dest="threads", default="8", help="Number of threads [8]")
parser.add_argument("-c", "--coverage", dest="coverage", default="30", help="Minimum coverage necessary for variant calling [30]")
parser.add_argument("-amb", "--ambiguity", dest="ambiguity", default="0", help="Define how to proceed with ambiguous genotypes. We offer 3 options: 0, no print; 1, print ambiguity codes; 2, randomly assign a nucleotide [0]")
parser.add_argument("-hapA", "--hapA", dest="hapA", action= "store", help="Haplotype A from reference genome (fasta). Warning: the name of each chromosome should be chr*_A in both haplotypes")
parser.add_argument("-hapB", "--hapB", dest="hapB", action= "store", help="Haplotype B from reference genome (fasta). Warning: the name of each chromosome should be chr*_B in both haplotypes")
parser.add_argument("-idA", "--TagName_hapA", dest="tagName_hapA", action="store", default="hapA", help="Tag to be added in the output files relative to hapA [default: hapA]")
parser.add_argument("-idB", "--TagName_hapB", dest="tagName_hapB", action="store", default="hapB", help="Tag to be added in the output files relative to hapB [default: hapB]")
parser.add_argument("-f1", "--fastq1", dest="fastq1", default="None", help="Illumina paired-end reads 1 (only required for read mapping)")
parser.add_argument("-f2", "--fastq2", dest="fastq2", default="None", help="Illumina paired-end reads 2 (only required for read mapping)")
parser.add_argument("-bA", "--bamA", dest="bamfileA", default="None", action="store", help="Bam file using haplotype A as reference")
parser.add_argument("-bB", "--bamB", dest="bamfileB", default="None", action="store", help="Bam file using haplotype B as reference")
parser.add_argument("-caller", "--variant_caller", dest="caller", default="GATK", help="Program for variant calling (GATK / bcftools / freebayes) [default: GATK]")
parser.add_argument("-coor", "--coordinatesTable", dest="coordinatesTable", default="None", help="Coordinates table. Required only if the two copies of each chromosome do not have a 1-to-1 position correspondence. Format: tab separated table with chr_hapA\tchr_HapB\tpositionA\tpositionB; ex: Ca22chrRA    Ca22chrRB    1    3")

args = parser.parse_args()


#############
#### TAG ####
#############

tagName_hapA = str(args.tagName_hapA)
tagName_hapB = str(args.tagName_hapB)
hapA = str(args.hapA)
hapB = str(args.hapB)
outDir = str(args.outDir)
fastq1 = str(args.fastq1)
fastq2 = str(args.fastq2)
numberthr = str(args.threads)
coverage = str(args.coverage)
coord = str(args.coordinatesTable)
amb = str(args.ambiguity)
bamA = str(args.bamfileA)
bamB = str(args.bamfileB)
caller = str(args.caller)

######################
#### All pipeline ####
######################
if fastq1 != "None" and fastq2 != "None" and bamA == "None" and bamB == "None":
	os.system(python_path+" " + haplotypo + "mapping.py -o " + outDir + " -thr " + numberthr + " -hapA " + hapA + " -hapB " + hapB + " -idA " + tagName_hapA  + " -idB " + tagName_hapB + " -f1 " + fastq1 + " -f2 " + fastq2)
	os.system(python_path+" " + haplotypo + "var_calling.py -o " + outDir + " -thr " + numberthr + " -hapA " + hapA + " -hapB " + hapB + " -idA " + tagName_hapA + " -idB " + tagName_hapB + " -bA " + outDir + "/" + tagName_hapA + ".mkd.bam " + " -bB " + outDir + "/" + tagName_hapB + ".mkd.bam " + " -c " + coverage + " -caller " + caller + " -p 2")
	if coord == "None":
		os.system(python_path+" " + haplotypo + "VCFcorr_alleles.py -A " + outDir + "/" + tagName_hapA + ".pass.snp.vcf" + " -B " + outDir + "/" + tagName_hapB + ".pass.snp.vcf" + " -fastaA " + hapA + " -fastaB " + hapB + " -cA " + outDir + "/" + tagName_hapA + ".corrected" + "_amb" + amb + ".vcf" + " -cB " + outDir + "/" + tagName_hapB + ".corrected" + "_amb" + amb + ".vcf" + " -amb " + amb) 
	else:
		os.system(python_path+" " + haplotypo + "VCFcorr_alleles.py -A " + outDir + "/" + tagName_hapA + ".pass.snp.vcf" + " -B " + outDir + "/" + tagName_hapB + ".pass.snp.vcf" + " -fastaA " + hapA + " -fastaB " + hapB + " -cA " + outDir + "/" + tagName_hapA + ".corrected" + "_amb" + amb + ".vcf" + " -cB " + outDir + "/" + tagName_hapB + ".corrected" + "_amb" + amb + ".vcf" + " -amb " + amb + " -coor " + coord) 
	os.system(python_path+" " + haplotypo + "haplomaker.py -o " + outDir + " -hapA " + hapA + " -hapB " + hapB + " -corrA " + outDir + "/" + tagName_hapA + ".corrected" + "_amb" + amb + ".vcf" + " -corrB " + outDir + "/" + tagName_hapB + ".corrected" + "_amb" + amb + ".vcf") 


######################################
#### Only from variant calling on ####
######################################
elif fastq1 == "None" and fastq2 == "None" and bamA != "None" and bamB != "None":
	os.system(python_path+" " + haplotypo + "var_calling.py -o " + outDir + " -thr " + numberthr + " -hapA " + hapA + " -hapB " + hapB + " -idA " + tagName_hapA + " -idB " + tagName_hapB + " -bA " + bamA + " -bB " + bamB + " -c " + coverage + " -caller " + caller + " -p 2")
	if coord == "None":
		os.system(python_path+" " + haplotypo + "VCFcorr_alleles.py -A " + outDir + "/" + tagName_hapA + ".pass.snp.vcf" + " -B " + outDir + "/" + tagName_hapB + ".pass.snp.vcf" + " -fastaA " + hapA + " -fastaB " + hapB + " -cA " + outDir + "/" + tagName_hapA + ".corrected" + "_amb" + amb + ".vcf" + " -cB " + outDir + "/" + tagName_hapB + ".corrected" + "_amb" + amb + ".vcf" + " -amb " + amb) 
	else:
		os.system(python_path+" " + haplotypo + "VCFcorr_alleles.py -A " + outDir + "/" + tagName_hapA + ".pass.snp.vcf" + " -B " + outDir + "/" + tagName_hapB + ".pass.snp.vcf" + " -fastaA " + hapA + " -fastaB " + hapB + " -cA " + outDir + "/" + tagName_hapA + ".corrected" + "_amb" + amb + ".vcf" + " -cB " + outDir + "/" + tagName_hapB + ".corrected" + "_amb" + amb + ".vcf" + " -amb " + amb + " -coor " + coord) 

	os.system(python_path+" " + haplotypo + "haplomaker.py -o " + outDir + " -hapA " + hapA + " -hapB " + hapB + " -corrA " + outDir + "/" + tagName_hapA + ".corrected" + "_amb" + amb + ".vcf" + " -corrB " + outDir + "/" + tagName_hapB + ".corrected" + "_amb" + amb + ".vcf") 
else:
	print("ERROR!! Please revise the fastq or bam parameters!")
