#!/usr/bin/env python

"""
Pipeline: Read mapping on phased or unphased genomes
12/02/2019
Veronica Mixao, Laia Carrete
CRG (Barcelona)
"""

import argparse
import os
import sys
import ConfigParser
from programs_config import program_path

#Setting the correct path for each program
bwa = program_path("bwa")
picard = program_path("picard")
samtools = program_path("samtools")

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
parser.add_argument("-f1", "--fastq1", dest="fastq1", required=True, help="Illumina paired-end reads 1")
parser.add_argument("-f2", "--fastq2", dest="fastq2", required=True, help="Illumina paired-end reads 2")

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

def mapping(reference, f1, f2, outdir, tag, threads):
	#Index the reference
	print "Indexing the reference..."
	cmd_indexFasta = bwa + " index " + reference
	print cmd_indexFasta
	os.system(cmd_indexFasta)
	
	#BWA-MEM
	print "\n", "Read mapping..."
	cmdbwa = bwa + ' mem -R "@RG\\tID:' + tag + '\\tSM:' + tag + '" -t ' + threads + " " + reference + " " + f1 + " " + f2 + " > " + outdir + "/" + tag + ".sam"
	print cmdbwa
	os.system(cmdbwa)
	print "BWA MEM done."
	
	#PICARD
	print "\n", "Sorting..."
	cmd_toBAM_sorted = picard + " SortSam " + "-I " + outdir + "/" + tag + ".sam " + "-O " + outdir + "/" + tag + ".bam" + " -SO coordinate" 
	print cmd_toBAM_sorted
	os.system(cmd_toBAM_sorted)
	
	print "\n", "Mark duplicates..."
	cmd_toBAM_mkd = picard + " MarkDuplicates " + "-I " + outdir + "/" + tag + ".bam " + "-O " + outdir + "/" + tag + ".mkd.bam" + " -M " + outdir + "/" + tag + ".mkd.metric.txt" 
	print cmd_toBAM_mkd
	os.system(cmd_toBAM_mkd)
	
	print "\n", "Indexing..."
	cmd_indexBAM = picard + " BuildBamIndex " + "-I " + outdir + "/" + tag + ".mkd.bam "
	print cmd_indexBAM
	os.system(cmd_indexBAM)
	
	print "Checking mapping quality..."
	cmd_qual = picard + " CollectAlignmentSummaryMetrics -R " + reference + " -I " + outdir + "/" + tag + ".mkd.bam " + " -O " + outdir + "/" + tag + ".mkd.stat.txt "
	print cmd_qual
	os.system(cmd_qual)
	
	print "\n", "Removing sam file to save space..."
	os.system("rm " + outdir + "/" + tag + ".bam") #this file is not needed anymore and takes too much space
	os.system("rm " + outdir + "/" + tag + ".sam") #this file is not needed anymore and takes too much space
	
	print "BAM file is done!"


#################
#### Running ####
#################

if ref == "None" and hapA != "None" and hapB != "None":	
	print "Read mapping will be performed on phased haplotypes..."
	print "*********************************************************************"
	print "Read mapping against haplotype A"
	mapping(hapA, fastq1, fastq2, outDir, tagName_hapA, numberthr)
	print "*********************************************************************"
	print "*********************************************************************"
	print "Read mapping against haplotype B"
	mapping(hapB, fastq1, fastq2, outDir, tagName_hapB, numberthr)
	print "*********************************************************************"
	
elif ref != "None" and hapA == "None" and hapB == "None":
	 print "Read mapping will be performed..."
	 mapping(ref, fastq1, fastq2, outDir, tagName, numberthr)
	
else:
	print "ERROR: Please revise the reference option you have given!"
