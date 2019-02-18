This repository has the Haplotypo.

This pipeline is still under development.

Haplotypo aims to do read mapping and variant calling on phased genomes, without loosing phasing information. In the end the corrected haplotypes are provided.

The different scripts represent:
1) mapping.py -> pipeline for read mapping with BWA-MEM
2) var_calling.py -> pipeline for SNP calling using one of the programs: GATK, bcftools or freebayes
3) VCFcorr_alleles.py -> script for to find the correspondence of SNPs between the two haplotypes and output a corrected VCF (only SNPs that actually belong to this haplotype)
4) haplomaker.py -> script to reconstruct the two haplotypes taking the corrected VCF
5) haplotypo.py -> pipeline to automatically run all the other pipelines
6) programs_config.py -> file to set the path of each program
