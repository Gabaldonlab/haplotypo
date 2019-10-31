#!/usr/bin/env python

"""
Pipeline: Variant calling on phased or unphased genomes
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
java = program_path("java")
samtools = program_path("samtools")
picard = program_path("picard")
gatk = program_path("gatk")
freebayes = program_path("freebayes")
bcftools = program_path("bcftools")
vcffilter = program_path("vcffilter")

##############
#### MAIN ####
##############

parser = argparse.ArgumentParser(description="Variant calling on phased or unphased genomes.")
parser.add_argument("-o", "--outdir", dest="outDir", action="store", help="Directory where the data will be stored")
parser.add_argument("-thr", "--threads", dest="threads", default="8", help="Number of threads [8]")
parser.add_argument("-hapA", "--hapA", dest="hapA", action= "store", default="None", help="Haplotype A from reference genome (fasta)")
parser.add_argument("-hapB", "--hapB", dest="hapB", action= "store", default="None", help="Haplotype B from reference genome (fasta)")
parser.add_argument("-ref", "--reference", dest="ref", action= "store", default="None", help="Unphased reference genome (fasta) - Only required if no phased haplotypes are provided")
parser.add_argument("-idA", "--TagName_hapA", dest="tagName_hapA", action="store", default="hapA", help="Tag to be added in the output files relative to hapA [default: hapA]")
parser.add_argument("-idB", "--TagName_hapB", dest="tagName_hapB", action="store", default="hapB", help="Tag to be added in the output files relative to hapB [default: hapB]")
parser.add_argument("-t", "--TagName", dest="tagName", action="store", default="strain", help="Tag to be added in the output files if using a non-phased reference [default: strain] - Only required if no phased haplotypes are provided")
parser.add_argument("-bA", "--bamA", dest="bamfileA", default="None", action="store", help="Bam file using haplotype A as reference")
parser.add_argument("-bB", "--bamB", dest="bamfileB", default="None", action="store", help="Bam file using haplotype B as reference")
parser.add_argument("-b", "--bam", dest="bamfile", default="None", action="store", help="Bam file using a non-phased reference - Only required if no phased haplotypes are provided")
parser.add_argument("-c", "--coverage", dest="coverage", default="30", help="Minimum coverage necessary for variant calling [30]")
parser.add_argument("-caller", "--variant_caller", dest="caller", default="GATK", help="Program for variant calling (GATK / bcftools / freebayes) [default: GATK]")
parser.add_argument("-p", "--ploidy", dest="ploidy", default="2", help="Ploidy [2]")
parser.add_argument("-pfile", "--ploidy_file", dest="ploidy_file", default="None", help="Ploidy file (only necessary for variant calling with bcftools and ploidy != 2")
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
numberthr = str(args.threads)
coverage = str(args.coverage)
caller = str(args.caller)
ploidy = str(args.ploidy)
bamA = str(args.bamfileA)
bamB = str(args.bamfileB)
bam = str(args.bamfile)
ploidy_file = str(args.ploidy_file)


#########################
#### Variant calling ####
#########################

#1. Indexing the reference
def index_reference(reference):
	print("Indexing the reference for variant calling...")
	cmd_index_ref = samtools + " faidx " + reference
	print(cmd_index_ref)
	os.system(cmd_index_ref)
	
	r = str(reference).split(".fa")
	print("\n", "Creating the reference dictionary...")
	cmd_dict = picard + " CreateSequenceDictionary" + " -R " + reference + " -O " + r[0] + ".dict"
	print(cmd_dict)
	os.system(cmd_dict)


#2. Variant calling	with GATK
def GATK(reference, outdir, tag, threads, cov, b, p):
	print("Running GATK...")
	cmd_gatk = java + " -jar " + gatk + " HaplotypeCaller -R " + reference + " -I " + b + " -O " + outdir + "/" + tag + ".vcf" + " -ploidy " + p + " --genotyping-mode DISCOVERY --standard-min-confidence-threshold-for-calling 30 -bamout " + outdir + "/" + tag + ".GATK.bam "
	print(cmd_gatk)
	os.system(cmd_gatk)

	print("\n", "Running GATK VariantFiltration...") 
	#Filtering with GATK recommendations for Haplotype Caller
	cmd_filter = java + " -jar " + gatk + " VariantFiltration -R " + reference + " -V " + outdir + "/" + tag + '.vcf -G-filter-name "heterozygous" -G-filter "isHet == 1" --filter-name "BadDepthofQualityFilter" -filter "DP <= ' + cov + ' || QD < 2.0 || MQ < 40.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" -O ' + outdir + "/" + tag + ".flt.vcf"
	print(cmd_filter)
	os.system(cmd_filter)

	#Separating PASS SNPs
	print("\n", "Separating SNPs...")
	cmd_sepSNP = java + " -jar " + gatk + " SelectVariants -R " + reference + " -V " + outdir + "/" + tag + ".flt.vcf" + " -select-type SNP -O " + outdir + "/" + tag + ".flt.snp.vcf"
	print(cmd_sepSNP)
	os.system(cmd_sepSNP)
	
	print("\n", "Getting PASS SNPs...")
	snp_file = open(outdir + "/" + tag + ".flt.snp.vcf", "r")
	snp = snp_file.readlines()
	
	passfile = open(outdir + "/" + tag + ".pass.snp.vcf" , "w+")
	
	for line in snp:
		l = line.split("\n")
		if "#" in line:
			print(l[0], file=passfile)
		else:
			if "PASS" in line:
				print(l[0], file=passfile)

	snp_file.close()
	passfile.close()

	print("\n", "Variant calling is done.")


#3. Variant calling with samtools -> We are missing separate SNPs and get the PASS file
def BCFTOOLS(reference, outdir, tag, threads, cov, b, p):
	print("Running Bcftools mpileup...")
	#A minimum mapping quality of 20 is being required (this is not default parameter but we want to reduce the chance of errors
	cmd_bcftools_mpileup = bcftools + " mpileup -a \"AD,DP\" -O b -f " + reference + " " + b  + " -o " + outdir + "/" + tag + ".mpileup.bcf" 
	print(cmd_bcftools_mpileup)
	os.system(cmd_bcftools_mpileup)
	
	print("\n", "Running Bcftools call...")
	#Bcftools only accepts ploidy different from 2 if we give a ploidy file
	if p == "None":
		cmd_bcftools = bcftools + " call -m -f GQ,GP -v -O v --threads " + threads + " -o " + outdir + "/" + tag + ".vcf " + outdir + "/" + tag + ".mpileup.bcf " 
	else:
		cmd_bcftools = bcftools + " call -m -f GQ,GP -v -O v --threads " + threads + " --ploidy-file " + p + " -o " + outdir + "/" + tag + ".vcf " + outdir + "/" + tag + ".mpileup.bcf " 
	print(cmd_bcftools)
	os.system(cmd_bcftools)
	
	#As there are no recommendations for bcftools, we decided to apply exclusively the filter for coverage. To apply harder filters please edit this command!
	#Only SNPs are PASS
	print("\n", "Filtering...")
	cmd_filter = bcftools + " filter -m x -e 'INFO/DP <= " + cov + "'" + " -O v --threads " + threads + ' -o ' + outdir + "/" + tag + ".flt.vcf " +  outdir + "/" + tag + ".vcf"
	print(cmd_filter)
	os.system(cmd_filter)
	
	#Separating PASS SNP file
	print("\n", "Getting PASS SNPs...") 
	cmd_get_snp = bcftools + " filter -m x -i \"TYPE = 'SNP'\" " + " -O v --threads " + threads + ' -o ' + outdir + "/" + tag + ".pass.snp.vcf " +  outdir + "/" + tag + ".flt.vcf"
	print(cmd_get_snp)
	os.system(cmd_get_snp)
	
	print("\n", "Variant calling is done.")
	
	
#4. Variant calling with freebayes
def FREEBAYES(reference, outdir, tag, threads, cov, b, p):
	print("Running freebayes...")
	#Minimum mapping quality set to 20 and minimum coverage is editable. Multi-nucleotides polimorphism calls are ignored, only INDELs and SNPs will be reported 
     #--haplotype-length -1 splits complex SNPs
	cmd_freebayes = freebayes + " -f " + reference + " -p " + p + " --min-coverage " + cov + " -b " + b + " --haplotype-length -1 -v " + outdir + "/" + tag + ".flt.vcf"
	print(cmd_freebayes)
	os.system(cmd_freebayes)
	
	print("\n", "Getting PASS SNPs...") 
	#Get only SNPs 
	cmd_pass_snps = vcffilter + ' -f "TYPE = snp & QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1 " ' + outdir + "/" + tag + ".flt.vcf > " + outdir + "/" + tag + ".pass.snp.vcf"     
	print(cmd_pass_snps)
	os.system(cmd_pass_snps) 
		
	print("\n", "Variant calling is done.")
	
	
#################
#### Running ####
#################
	
if ref == "None" and hapA != "None" and hapB != "None":
	print("Variant calling will be performed on phased haplotypes...")
	if bam != "None":
		print("bam argument was not expected!")
	else:
		print("Variant calling will be performed...")	
		if caller == "GATK":
			index_reference(hapA)
			index_reference(hapB)
			print("*********************************************************************")
			print("Variant calling against haplotype A")
			GATK(hapA, outDir, tagName_hapA, numberthr, coverage, bamA, ploidy)
			print("*********************************************************************")
			print("*********************************************************************")
			print("Variant calling against haplotype B")
			GATK(hapB, outDir, tagName_hapB, numberthr, coverage, bamB, ploidy)
			print("*********************************************************************")
		elif caller == "bcftools":
			print("*********************************************************************")
			print("Variant calling against haplotype A")
			BCFTOOLS(hapA, outDir, tagName_hapA, numberthr, coverage, bamA, ploidy_file)
			print("*********************************************************************")
			print("*********************************************************************")
			print("Variant calling against haplotype B")
			BCFTOOLS(hapB, outDir, tagName_hapB, numberthr, coverage, bamB, ploidy_file)
			print("*********************************************************************")
		elif caller == "freebayes":
			print("*********************************************************************")
			print("Variant calling against haplotype A")
			FREEBAYES(hapA, outDir, tagName_hapA, numberthr, coverage, bamA, ploidy)
			print("*********************************************************************")
			print("*********************************************************************")
			print("Variant calling against haplotype B")
			FREEBAYES(hapB, outDir, tagName_hapB, numberthr, coverage, bamB, ploidy)
			print("*********************************************************************")
		else:
			print("Please provide a valid caller!")
elif ref != "None" and hapA == "None" and hapB == "None":
	if bamA != "None" or bamB != "None":
		print("bamA and bamB arguments are not expected!")
	else:
		print("Variant calling will be performed...")
		if caller == "GATK":
			index_reference(ref)
			GATK(ref, outDir, tagName, numberthr, coverage, bam, ploidy)
		elif caller == "bcftools":
			if ploidy != "2" and ploidy_file == "None":
				print("Please provide a ploidy file if ploidy != 2!")
			else:
				BCFTOOLS(ref, outDir, tagName, numberthr, coverage, bam, ploidy_file)
		elif caller == "freebayes":
			FREEBAYES(ref, outDir, tagName, numberthr, coverage, bam, ploidy)
		else:
			print("Please provide a valid caller!")
else:
	print("ERROR: Please revise the reference option you have given!")
