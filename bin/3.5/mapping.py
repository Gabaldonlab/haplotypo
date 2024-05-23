#!/usr/bin/env python

"""
Pipeline: Read mapping on phased or unphased genomes
23/05/2024
Veronica Mixao, Laia Carrete, Alvaro Redondo
CRG (Barcelona) - BSC
"""

import argparse
import os
import sys
import configparser
from programs_config import program_path

#Setting the correct path for each program
bwa = program_path("bwa")
picard = program_path("picard")

##############
#### MAIN ####
##############

parser = argparse.ArgumentParser(description="Read mapping on phased or unphased genomes.")
parser.add_argument("-o", "--outdir", dest="outDir", action="store", help="Directory where the data will be stored")
parser.add_argument("-thr", "--threads", dest="threads", default="8", help="Number of threads [8]")
parser.add_argument("-hapA", "--hapA", dest="hapA", action= "store", default="None", help="Haplotype A from reference genome (fasta)")
parser.add_argument("-hapB", "--hapB", dest="hapB", action= "store", default="None", help="Haplotype B from reference genome (fasta)")
parser.add_argument("-ref", "--reference", dest="ref", action= "store", default="None", help="Unphased reference genome (fasta) - Only required if no phased haplotypes are provided")
parser.add_argument("-idA", "--TagName_hapA", dest="tagName_hapA", action="store", default="hapA", help="Tag to be added in the output files relative to hapA [default: hapA]")
parser.add_argument("-idB", "--TagName_hapB", dest="tagName_hapB", action="store", default="hapB", help="Tag to be added in the output files relative to hapB [default: hapB]")
parser.add_argument("-t", "--TagName", dest="tagName", action="store", default="strain", help="Tag to be added in the output files if using a non-phased reference [default: strain] - Only required if no phased haplotypes are provided")
parser.add_argument("-f1", "--fastq1", dest="fastq1", required=True, help="Illumina paired-end reads 1 / Single-end reads")
parser.add_argument("-f2", "--fastq2", dest="fastq2", default="None", help="Illumina paired-end reads 2")
parser.add_argument("-s", "--single_end", dest="single_end", action="store_true", help="Flag for single-end reads")

args = parser.parse_args()


#############
#### TAG ####
#############

tagName_hapA = str(args.tagName_hapA)
tagName_hapB = str(args.tagName_hapB)
tagName = str(args.tagName)
hapA = str(args.hapA)
hapB = str(args.hapB)
ref = str(args.ref)
outDir = str(args.outDir)
fastq1 = str(args.fastq1)
fastq2 = str(args.fastq2)
numberthr = str(args.threads)


######################
#### Read Mapping ####
######################

def mapping(reference, single_end, outdir, tag, threads, f1, f2 = "None"):
	#Index the reference
	print("Indexing the reference...")
	cmd_indexFasta = f"{bwa} index {reference}"
	print(cmd_indexFasta)
	os.system(cmd_indexFasta)
	
	#BWA-MEM
	print("\n", "Read mapping...")
	if single_end:
		print("Mapping single-end reads")
		cmdbwa = (
    		f"{bwa} mem "
			f"-R \"@RG\\tID:{tag}\\tSM:{tag}\" -t {threads} "
    		f"{reference} {f1} "
			f"> {outdir}/{tag}.sam"
		)
	else:
		print("Mapping paired-end reads")
		cmdbwa =(
    		f"{bwa} mem "
			f"-R \"@RG\\tID:{tag}\\tSM:{tag}\" -t {threads} "
    		f"{reference} {f1} {f2} "
    		f"> {outdir}/{tag}.sam"
		)
	print(cmdbwa)
	os.system(cmdbwa)
	print("BWA MEM done.")
	
	#PICARD
	print("\n", "Sorting...")
	cmd_toBAM_sorted = (
		f"{picard} SortSam "
		f"-I {outdir}/{tag}.sam "
		f"-O {outdir}/{tag}.bam "
		f"-SO coordinate "
	)
	print(cmd_toBAM_sorted)
	os.system(cmd_toBAM_sorted)
	
	print("\n", "Mark duplicates...")
	cmd_toBAM_mkd = (
		f"{picard} MarkDuplicates "
		f"-I {outdir}/{tag}.bam "
		f"-O {outdir}/{tag}.mkd.bam "
		f"-M {outdir}/{tag}.mkd.metric.txt"
	)
	print(cmd_toBAM_mkd)
	os.system(cmd_toBAM_mkd)
	
	print("\n", "Indexing...")
	cmd_indexBAM = (
		f"{picard} BuildBamIndex "
		f"-I {outdir}/{tag}.mkd.bam"
	)
	print(cmd_indexBAM)
	os.system(cmd_indexBAM)
	
	print("Checking mapping quality...")
	cmd_qual = (
		f"{picard} CollectAlignmentSummaryMetrics "
		f"-R {reference} "
		f"-I {outdir}/{tag}.mkd.bam "
		f"-O {outdir}/{tag}.mkd.stat.txt"
	)
	print(cmd_qual)
	os.system(cmd_qual)
	
	print("\n", "Removing sam file to save space...")
	os.system("rm " + outdir + "/" + tag + ".bam") #this file is not needed anymore and takes too much space
	os.system("rm " + outdir + "/" + tag + ".sam") #this file is not needed anymore and takes too much space
	
	print("BAM file is done!")


#################
#### Running ####
#################

if not args.single_end and fastq2 == "None":
	parser.error("ERROR: Reverse reads not found, use the -s flag for single-end data")

if ref == "None" and hapA != "None" and hapB != "None":	
	print("Read mapping will be performed on phased haplotypes...")
	print("*********************************************************************")
	print("Read mapping against haplotype A")
	mapping(hapA, args.single_end, outDir, tagName_hapA, numberthr, fastq1, fastq2)
	print("*********************************************************************")
	print("*********************************************************************")
	print("Read mapping against haplotype B")
	mapping(hapB, args.single_end, outDir, tagName_hapB, numberthr, fastq1, fastq2)
	print("*********************************************************************")
	
elif ref != "None" and hapA == "None" and hapB == "None":
	print("Read mapping will be performed...")
	mapping(ref, args.single_end, outDir, tagName, numberthr, fastq1, fastq2)

	
else:
	print("ERROR: Please revise the reference option you have given!")
