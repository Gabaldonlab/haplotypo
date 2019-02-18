#!/usr/bin/env python

"""
Pipeline: Reconstructing the phased haplotypes for the strain in analysis
12/02/2019
Veronica Mixao
CRG (Barcelona)

python haplomaker.py -o output_directory -hapA hapA.fa -hapB hapB.fa -corrA strain_hapA.corrected_amb0.vcf -corrB strain_hapB.corrected_amb0.vcf
"""

import argparse


##############
#### MAIN ####
##############

parser = argparse.ArgumentParser(description="Reconstruction of phased haplotypes based on VCFcorr_DEF.py output.")
parser.add_argument("-o", "--outdir", dest="outDir", action="store", help="Directory where the data will be stored")
parser.add_argument("-hapA", "--hapA", dest="hapA", action= "store", default="None", help="Haplotype A from reference genome (fasta)")
parser.add_argument("-hapB", "--hapB", dest="hapB", action= "store", default="None", help="Haplotype B from reference genome (fasta)")
parser.add_argument("-corrA", "--corrected_A", dest="corrected_A", action="store", help="VCF file with the corrected variants for haplotype A")
parser.add_argument("-corrB", "--corrected_B", dest="corrected_B", action="store", help="VCF file with the corrected variants for haplotype B")

args = parser.parse_args()


#############
#### TAG ####
#############

hapA = str(args.hapA)
hapB = str(args.hapB)
outDir = str(args.outDir)
corrA = str(args.corrected_A)
corrB = str(args.corrected_B)


#######################################
#### Reconstructing the haplotypes ####
#######################################

def haploMaker(fasta,vcf):
	f_open = open(fasta, "r")
	vcf_open = open(vcf, "r")
	
	f = f_open.readlines()
	var = vcf_open.readlines()
	
	ref = {}

	#Getting the reference sequence (ref[chrm] = seq)
	for line in f:
		if ">" in line:
			seq = ""
			lin = line.split("\n")
			l = lin[0].split(">")
			chrm = l[1]
		else:
			l = line.split("\n")
			seq += l[0]
		ref[chrm] = seq
	
	#Getting information from the new VCF file (alt_nucl[chrm][pos] = ref,alt)
	alt_nucl = {}
	pos_order_chromosome = {}

	for line in var:
		if "#" not in line:
			l = line.split("\t")
			chrm = l[0]
			pos = l[1]
			initial = l[3]
			alt = l[4]
			
			infor = initial,alt
			
			if chrm not in alt_nucl.keys():
				alt_nucl[chrm] = {}
				pos_order_chromosome[chrm] = []
			alt_nucl[chrm][pos] = infor
			pos_order_chromosome[chrm].append(pos)
	
	#Changing the haplotype reference according to SNP information
	fasta_name = vcf.split(".")
	
	if "/" in fasta_name[0]:
		name = fasta_name[0].split("/")
		tag = name[len(name)-1]
	else:
		tag = fasta_name[0]
		
	new_fasta = open(outDir + "/" + tag + ".alternative.fasta", "w+")
	
	for chromosome in ref.keys():
		seq = "" #new sequence
		last_position = 0
		if chromosome in pos_order_chromosome.keys():
			alternative_positions = pos_order_chromosome[chromosome]
			for position in alternative_positions:
				seq += ref[chromosome][int(last_position):int(position)-1]
				initial,alt = alt_nucl[chromosome][position]
				if initial == ref[chromosome][int(position)-1]:
					if len(alt) == 1: #if only one alternative, no ambiguity code needed
						code = alt
					else:
						if len(alt) >= 2: #if more than one alternative, ambiguity code needed
							print "WARNING: More than one alternative possibility in the corrected VCF!", chromosome, position, initial, alt
					seq += code #add the new nucleotide to the reference
				else:
					print "WARNING!!!! The nucleotide in VCF does not correspond to nucleotide in the reference" #print a warning because the reference nucleotide has to match the reference in the VCF for the same position
					print chromosome, position, initial, ref[chromosome][int(position)-1], alt
				last_position = int(position)
			if last_position < len(ref[chromosome]):
				end_chromosome = ref[chromosome][int(last_position):] #after the last alternative position we need to continue adding the nucleotides in case it does not represent the last position of the chromosome
				seq += end_chromosome
		else:
			seq += ref[chromosome]
		print >>new_fasta, ">" + chromosome + "_alternative" + "\n" + seq

	f_open.close()
	vcf_open.close()
	new_fasta.close()

print "Generating the corrected haplotype A..."
haploMaker(hapA,corrA)
print "Generating the corrected haplotype B..."
haploMaker(hapB,corrB)

print "\n", "Haplomaker is done."
