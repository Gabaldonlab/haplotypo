# Haplotypo
![Docker Pulls](https://img.shields.io/docker/pulls/cgenomics/haplotypo)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://hub.docker.com/repository/docker/cgenomics/haplotypo)

This repository contains the Haplotypo pipeline.

## Description

HaploTypo is a pipeline suited to map variants into haplotypes in genetic variation analyses. After mapping and variant calling on a phased reference genome, HaploTypo infers the haplotype correspondence for each heterozygous variant. It also generates two independent FASTA files for each reconstructed haplotype.

## Scripts contained in this repository:
1) mapping.py -> pipeline for read mapping with BWA-MEM
2) var_calling.py -> pipeline for SNP calling using one of the programs: GATK, bcftools or freebayes
3) VCFcorr_alleles.py -> script for to find the correspondence of SNPs between the two haplotypes and output a corrected VCF (only SNPs that actually belong to this haplotype)
4) haplomaker.py -> script to reconstruct the two haplotypes taking the corrected VCF
5) haplotypo.py -> pipeline to automatically run all the other pipelines
6) programs_config.py -> file to set the path of each program
