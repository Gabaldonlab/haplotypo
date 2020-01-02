# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 12:45:14 2018

@author: cpegueroles
"""
#sudo pip install pyvcf

import argparse, pyvcf, re, random, sys,os
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Pipeline to obtain corrected VCF files for haplotypes A and B from their filtered VCFs. It requires pyVCF. WARNING: chromosome Id must be /(chr.*?)_.*/; ex: chrRA")

parser.add_argument("-A", "--VCF_A", dest="VCF_A", required=True, help="VCF with PASS SNPs obtained after running a variant caller using HapA as reference")
parser.add_argument("-B", "--VCF_B", dest="VCF_B", required=True, help="VCF with PASS SNPs obtained after running a variant caller using HapB as reference")
parser.add_argument("-fastaA", "--fasta_A", dest="fasta_A", required=True, help="Fasta with HapA used as reference for variant calling")
parser.add_argument("-fastaB", "--fasta_B", dest="fasta_B", required=True, help="Fasta with HapB used as reference for variant calling")
parser.add_argument("-cA", "--VCF_corrA", dest="VCF_corrA", required=True, help="Out file for the corrected VCF A")
parser.add_argument("-cB", "--VCF_corrB", dest="VCF_corrB", required=True, help="Out file for the corrected VCF B")
parser.add_argument("-amb", "--ambiguity", dest="ambiguity", required=True, help="Define how to proceed with ambiguous genotypes-> 0: no print; 1: print ambiguity codes; 2: randomly assign a ncl")
parser.add_argument("-coor", "--coordinatesTable", dest="coordinatesTable", default="None", help="Coordinates table. Required only if the two copies of each chromosome have different lengths. Format: tab separated table with chr_hapA\tchr_HapB\tpositionA\tpositionB; ex: chrRA	chrRB	1	1")


args = parser.parse_args()
VCF_rawA = str(args.VCF_A)
VCF_rawB = str(args.VCF_B)
fasta_A = str(args.fasta_A)
fasta_B = str(args.fasta_B)
VCF_corrA = str(args.VCF_corrA)
VCF_corrB = str(args.VCF_corrB)
ambiguity = str(args.ambiguity)

if ambiguity != '0' and ambiguity != '1' and ambiguity != '2':
    sys.exit("\tWARNING! Please check your --amb option, should be 0, 1 or 2")
        
########## parse VCF header ############
def parse_vcf_header(VCF_raw,VCF_corr,ambiguity): #parse header and print the info in the output file
    VCF_rawFH = open(VCF_raw, 'r')
    for line in VCF_rawFH:
        line = line.rstrip()
        if line.startswith('##fileformat') and ambiguity != '1':
            VCF_corr.write('%s\n' % line)
        if line.startswith('##fileformat') and ambiguity == '1': 
            VCF_corr.write('##fileformat=VCFv4.3\n')
        if line.startswith('##FILTER'):
            VCF_corr.write('%s\n' % line)
        if line.startswith('##FORMAT=<ID=AD'):
            VCF_corr.write('%s\n' % line)
        if line.startswith('##FORMAT=<ID=DP'):
            VCF_corr.write('%s\n' % line)
        if line.startswith('##FORMAT=<ID=GT'):
            VCF_corr.write('%s\n' % line)
        if line.startswith('##contig'):
            VCF_corr.write('%s\n' % line)
        if line.startswith('#CHROM'):
            if ambiguity == '1':
                VCF_corr.write('##ALT=<ID=Y,Description="IUPAC code Y = T/C">\n##ALT=<ID=W,Description="IUPAC code W = T/A">\n##ALT=<ID=R,Description="IUPAC code R = A/G">\n##ALT=<ID=S,Description="IUPAC code S = C/G">\n##ALT=<ID=K,Description="IUPAC code K = T/G">\n##ALT=<ID=M,Description="IUPAC code M = A/C">\n')   
            VCF_corr.write('%s\n' % line)            
            
########## VCF read if coordinates table is NOT provided #################
def read_vcf(VCF_raw):
    vcf_dict={}
    vcf_reader = pyvcf.Reader(open(VCF_raw, 'r'))
    for record in vcf_reader:
        #print record.CHROM,record.POS,record.ID,record.REF,record.ALT,record.QUAL,record.FILTER,record.INFO,record.samples
        m = re.search('(chr.*?)_.*', record.CHROM)
        mychr=m.group(1)
        #print mychr,record.POS
        if not mychr:
            print("WARNING! Chromosome Id must be /(chr.*?)_.*/; ex: chr1_A; chr1_B; please check! ;)")
        key= "%s_%s" % (mychr,record.POS)
        val= record
        if not key in vcf_dict:
            vcf_dict[key]=val
        else:
            print("WARNING: %s  and %s appear more than once" % (record.CHROM,record.POS))
    return vcf_dict
    
########## VCF read if coordinates table is provided #################
def read_vcfDiffCoord(VCF_raw,haplotype, coord): ####new
    vcf_dict={}
    vcf_reader = pyvcf.Reader(open(VCF_raw, 'r'))
    for record in vcf_reader:
        m = re.search('(chr.*?)_.*', record.CHROM)
        mychr=m.group(1)
        newPos=int()
        found=0
        if haplotype=='A':
            for x in coord:
                if x[0]==mychr and int(x[2])==int(record.POS):
                    newchr=x[1]
                    newPos=x[3]
                    found=1
            key= "%s_%s" % (newchr,newPos)
        if found ==0:
            #for positions in hapA with no correspondance in the coordinates table it recodes to position + length of the chr
            #the reason is that this certain position in A can exist in B depite not being its correspondance
            print ("\tWARNING: please check your coordinated file, position %s from chromosome %s is missing!" % (record.POS,m.group(0)))
            for fasta in SeqIO.parse(fasta_A, "fasta"):
                if fasta.id == '%s_A' %m.group(1):
                   mylen=len(fasta.seq)
            key= "%s_%s" % (mychr,(int(record.POS)+int(mylen)))
        val= record
        if not key in vcf_dict:
            vcf_dict[key]=val
        else:
            print("WARNING: %s  and %s appear more than once" % (record.CHROM,record.POS))
    return vcf_dict

########## save Genotype in dict #################   
def save_GT(GT1_corr, key1, val1,vcfcorr):
    key=key1
    val = {'chr':[], 'position':[], 'reference':[], 'alternative':[],'GT':[], 'AD':[], 'DP':[]}
    if not key in vcfcorr:
        vcfcorr[key] =val
    vcfcorr[key]['chr'].append(val1.CHROM)
    vcfcorr[key]['position'].append(val1.POS)
    vcfcorr[key]['reference'].append(val1.REF)
    vcfcorr[key]['alternative'].append(val1.ALT[0])
    vcfcorr[key]['GT'].append(GT1_corr)
    vcfcorr[key]['AD'].append(val1.samples[0]['AD'])                      
    vcfcorr[key]['DP'].append(val1.samples[0]['DP'])
    return (vcfcorr)
    
########## ambiguity code #################   
def ambGenotypes(alt):
    ambDict = { "Y": ["T", "C"], "R": ["A", "G"], "W": ["T", "A"], "S": ["C", "G"] , "K": ["T", "G"], "M": ["C", "A"] }
    for key, val in ambDict.items():
        nalt=[str(alt[0]),str(alt[1])]
        if sorted(nalt) == sorted(val):
            return (key)    
            
########## random genotypes ################# 
def randomGenotypes(alt, randomDict, key1):
    key=key1
    val=''
    if not key in randomDict:
        randomDict[key]=val
        randomDict[key] = alt[random.randint(0, 1)] #print the alternative or reference alleles randomly
        newAlt=str(randomDict[key])
    else:
        gt=[x for x in alt if str(x) != str(randomDict[key])]
        newAlt=str(gt[0])
    return [newAlt, randomDict]  
    
########## ambiguity code for 0/1-0/1 #################   
def ambGenotypes01(ref,alt):
    ambDict = { "Y": ["T", "C"], "R": ["A", "G"], "W": ["T", "A"], "S": ["C", "G"] , "K": ["T", "G"], "M": ["C", "A"] }
    for key, val in ambDict.items():
        nalt=[str(ref),str(alt[0])]     
        if sorted(nalt) == sorted(val):
            return (key)    
            
########## random genotypes for 0/1-0/1 ################# 
def randomGenotypes01(ref,alt, randomDict, key1):
    key=key1
    val=''
    nalt=[str(ref),str(alt[0])]
    if not key in randomDict:
        randomDict[key]=val
        randomDict[key] = nalt[random.randint(0, 1)] #print the alternative or reference alleles randomly
        newAlt=str(randomDict[key])
    else:
        gt=[x for x in nalt if str(x) != str(randomDict[key])]
        newAlt=str(gt[0])
    return [newAlt, randomDict]  
    
########## get refence ncl when no SNP is annotated ################# 
def getRefNcl(key1,haplotype):
    chrPos=key1.split('_')
    if haplotype == 'A':
        for seq_record in SeqIO.parse(fasta_B, "fasta"):
            m = re.search('chr.*?_(.*)', seq_record.id)
            tag=m.group(1)
            if seq_record.id == '%s_%s' % (chrPos[0],tag): #get the chromosome from the haplotype that do not have SNPs annotated
                return seq_record.seq[int(chrPos[1])-1]
    if haplotype == 'B':
        for seq_record in SeqIO.parse(fasta_A, "fasta"):
            m = re.search('chr.*?_(.*)', seq_record.id)
            tag=m.group(1)
            if seq_record.id == '%s_%s' % (chrPos[0],tag):#'%s_A'% chrPos[0]: #get the chromosome from the haplotype that do not have SNPs annotated
                return seq_record.seq[int(chrPos[1])-1]    
    
########## VCF correct #################   
def VCF_corr(vcf1_dict,vcf2_dict,randomDict,haplotype, coord):
    vcfcorr={}
    unsolvedGT=0
    unsolved=[]; newUnphasedVariants=[]
    for key1, val1 in vcf1_dict.items():
        GT1_corr=''; found=0
        if key1 in vcf2_dict:
            found01_01=0; found11_11=0; found01_12=0; found12_01=0; found12_12=0
            if val1.samples[0]['GT']=='0/1' and vcf2_dict[key1].samples[0]['GT']=='0/1': #unphased 
                if str(val1.REF) ==str(vcf2_dict[key1].REF) and str(val1.ALT[0])==str(vcf2_dict[key1].ALT[0]):
                    found01_01=1
                    newUnphasedVariants.append([val1.CHROM, val1.POS])
                    if ambiguity == '0':
                        GT1_corr=0                       
                    if ambiguity == '1':
                        GT1_corr=1
                        save_GT(GT1_corr, key1, val1,vcfcorr)
                        vcfcorr[key1]['alternative']=ambGenotypes01(val1.REF,val1.ALT)#to store an ambiguity code instead of two possibilities
                    if ambiguity == '2':
                        res=randomGenotypes01(val1.REF,val1.ALT, randomDict, key1)# to randomly assign one of the two possible ncl                            
                        if val1.REF!=res[0]:
                            GT1_corr=1
                        else:
                            GT1_corr=0
                        save_GT(GT1_corr, key1, val1,vcfcorr)
                        vcfcorr[key1]['alternative']=res[0]
                        randomDict=res[1]
                if str(val1.REF) ==str(vcf2_dict[key1].ALT[0]) and str(vcf2_dict[key1].REF)==str(val1.ALT[0]): #solved 
                    found01_01=1
                    GT1_corr=0
                    save_GT(GT1_corr, key1, val1,vcfcorr)
                if found01_01==0 : #unsolved
                    unsolvedGT+=1
                    unsolved.append([val1.CHROM, val1.POS])
                found=1
            if val1.samples[0]['GT']=='1/1' and vcf2_dict[key1].samples[0]['GT']=='1/1':
                if str(val1.ALT[0])==str(vcf2_dict[key1].ALT[0]): #solved
                    found11_11=1
                    GT1_corr=1
                    #print 'found 1/1-1/1'
                    save_GT(GT1_corr, key1, val1,vcfcorr)
                if found11_11==0: #unsolved
                    unsolvedGT+=1
                    unsolved.append([val1.CHROM, val1.POS])
                found=1
            if val1.samples[0]['GT']=='0/1' and vcf2_dict[key1].samples[0]['GT']=='1/2':  
                if str(val1.REF) in str(vcf2_dict[key1].ALT) and str(val1.ALT[0]) in str(vcf2_dict[key1].ALT):
                    found01_12=1
                    GT1_corr=0
                    save_GT(GT1_corr, key1, val1,vcfcorr)
                if found01_12==0: #unsolved
                    unsolvedGT+=1
                    unsolved.append([val1.CHROM, val1.POS])
                found=1
            if val1.samples[0]['GT']=='1/2' and vcf2_dict[key1].samples[0]['GT']=='0/1':
                if str(vcf2_dict[key1].REF) in str(val1.ALT) and str(vcf2_dict[key1].ALT[0]) in str(val1.ALT): #now the order is reversed compared to 0/1-1/2 
                    found12_01=1
                    GT1_corr=1
                    allele1=[x for x in val1.ALT if x in vcf2_dict[key1].ALT]
                    allele2=[x for x in val1.ALT if x not in vcf2_dict[key1].ALT]
                    val1.ALT=[allele1[0],allele2[0]] #reorder the alleles since in save_GT prints the first allele
                    save_GT(GT1_corr, key1, val1,vcfcorr)
                if found12_01==0: #unsolved
                    unsolvedGT+=1
                    unsolved.append([val1.CHROM, val1.POS])
                found=1
            if val1.samples[0]['GT']=='1/2' and vcf2_dict[key1].samples[0]['GT']=='1/2':
                if sorted([str(val1.ALT[0]),str(val1.ALT[1])]) == sorted([str(vcf2_dict[key1].ALT[0]),str(vcf2_dict[key1].ALT[1])]):#unphased
                    found12_12=1
                    newUnphasedVariants.append([val1.CHROM, val1.POS])                        
                    if ambiguity == '0':
                            GT1_corr=0                       
                    if ambiguity == '1':
                        GT1_corr=1
                        save_GT(GT1_corr, key1, val1,vcfcorr)
                        vcfcorr[key1]['alternative']=ambGenotypes(val1.ALT)#to store an ambiguity code instead of two possibilities
                    if ambiguity == '2':
                        res=randomGenotypes(val1.ALT, randomDict, key1)# to randomly assign one of the two possible ncl
                        if val1.REF!=res[0]:
                            GT1_corr=1
                        else:
                            GT1_corr=0
                        save_GT(GT1_corr, key1, val1,vcfcorr)
                        vcfcorr[key1]['alternative']=res[0]
                        randomDict=res[1]                          
                if found12_12==0 : #unsolved
                    unsolvedGT+=1
                    unsolved.append([val1.CHROM, val1.POS])
                found=1
            if val1.samples[0]['GT']=='0/1' and vcf2_dict[key1].samples[0]['GT']=='1/1':
                unsolvedGT+=1
                unsolved.append([val1.CHROM, val1.POS])
                found=1
            if val1.samples[0]['GT']=='1/1' and vcf2_dict[key1].samples[0]['GT']=='0/1':
                unsolvedGT+=1
                unsolved.append([val1.CHROM, val1.POS])
                found=1
            if val1.samples[0]['GT']=='1/1' and vcf2_dict[key1].samples[0]['GT']=='1/2':
                unsolvedGT+=1
                unsolved.append([val1.CHROM, val1.POS])
                found=1
            if val1.samples[0]['GT']=='1/2' and vcf2_dict[key1].samples[0]['GT']=='1/1':
                unsolvedGT+=1
                unsolved.append([val1.CHROM, val1.POS])
                found=1       
                
        if val1.samples[0]['GT']=='1/1' and found==0:
            found11_NoSNP=0
            key1old=''; positionFoundinCoorTablB=0; positionCheckA=0
            m = re.search('(chr.*?)_(.*)', key1) 
            if ctable==0:
                if getRefNcl(key1, haplotype)==val1.ALT[0]: #solved
                    found11_NoSNP=1
                    GT1_corr=1
                    save_GT(GT1_corr, key1, val1,vcfcorr)
            if ctable==1 and haplotype=='B':  #to check which position get from haplotype A 
                key1old=key1
                positionFoundinCoorTablB=1
                for x in coord:
                    if m.group(1)==x[1] and m.group(2)==x[3]:
                        key1='%s_%s' %(x[0],x[2]) 
                        positionFoundinCoorTablB=2
            for record in SeqIO.parse(fasta_A, "fasta"):
                if record.id == '%s_A' %m.group(1):
                   mylen=record.seq
            if ctable==1 and haplotype=='A': 
                key1old=key1
                positionCheckA=1
                if int(m.group(2)) <= int(len(mylen)): 
                    positionCheckA=2 
            if positionCheckA==2 or positionFoundinCoorTablB==2:
                if getRefNcl(key1, haplotype)==val1.ALT[0]: #solved
                    found11_NoSNP=1
                    GT1_corr=1
                    save_GT(GT1_corr, key1old, val1,vcfcorr)
            if positionCheckA==1 or positionFoundinCoorTablB==1:
                    found11_NoSNP=1
                    GT1_corr=1
                    save_GT(GT1_corr, key1old, val1,vcfcorr)
            if found11_NoSNP==0 :
                unsolvedGT+=1
                unsolved.append([val1.CHROM, val1.POS])               
        if val1.samples[0]['GT']=='0/1' and found==0: #val1.samples[0]['GT']=='0/1' and no SNP in vcf2_dict[key1]; biologically not possible, is the results of the filtering
            unsolvedGT+=1
            unsolved.append([val1.CHROM, val1.POS])
        if val1.samples[0]['GT']=='1/2' and found==0: #val1.samples[0]['GT']=='1/2' and no SNP in vcf2_dict[key1]; biologically not possible, is the results of the filtering
            unsolvedGT+=1
            unsolved.append([val1.CHROM, val1.POS])  
            
    return [vcfcorr,randomDict,unsolved,newUnphasedVariants]  #I need to return randomDict
##########################################

######## to read the raw VCF files ##################################
coord=[]
if args.coordinatesTable is 'None':
    print("\n\tcoordinates table is NOT provided")
    vcfA_dict= read_vcf(VCF_rawA)
    vcfB_dict= read_vcf(VCF_rawB)
    ctable=0 

else:
    print("\n\tcoordinates table IS provided")
    coordinatesTable= str(args.coordinatesTable)
    for line in open(coordinatesTable, 'r'):
        line = line.rstrip().split('\t')
        m = re.search('(chr.*?)_.*', line[0])
        chrA=m.group(1)
        m = re.search('(chr.*?)_.*', line[1]) 
        chrB=m.group(1)
        coord.append([chrA,chrB,line[2],line[3]]) #to save the coordinates in a list of lists 

    vcfA_dict= read_vcfDiffCoord(VCF_rawA,'A',coord) 
    vcfB_dict= read_vcf(VCF_rawB)
    ctable=1 


randomDict={} #to initialize randomDict

######## to obtain and print corrected VCF from haplotype A ##################################
print('\n\thaplotype A')
res= VCF_corr(vcfA_dict,vcfB_dict,randomDict,'A', coord) #obtain the corrected VCF; WARNING: the order of the dic is important!! ####new
vcfA_corr_dict=res[0]
randomDict=res[1]
unsolvedA=res[2]
unphasedA=res[3]

tag=os.path.abspath(args.VCF_A)
m=re.search('(.*)\.pass\.snp\.vcf', tag)
tag=m.group(1)
unsolvedAout=open('%s.unsolved_amb%s.bed' % (tag, ambiguity), 'w')
unphasedAout=open('%s.unphasedVariants_amb%s.bed' % (tag, ambiguity), 'w')


for x in sorted(unsolvedA, key=lambda x: (x[0], int(x[1]))): #to print SORTED chr and position of unsolved cases (see VCF_corr function)
    unsolvedAout.write('%s\t%s\t%s\n' % (x[0],int(x[1])-1,x[1]))

for x in sorted(unphasedA, key=lambda x: (x[0], int(x[1]))): #to print SORTED chr and position of unphased cases (see VCF_corr function)
    unphasedAout.write('%s\t%s\t%s\n' % (x[0],int(x[1])-1,x[1]))
    
VCF_corrA= open(VCF_corrA, 'w')
parse_vcf_header(VCF_rawA,VCF_corrA,ambiguity) #parse and write header in VCF corr

### to print the sorted vcf file:
myList=[]
for key, val in vcfA_corr_dict.items():
    if val['GT'][0] ==1:#to append only positions with GT=1
        myList.append([val['chr'][0], val['position'][0],'.',val['reference'][0], val['alternative'][0],'.','.','.','GT:AD:DP\t%s:%s,%s:%s' % (val['GT'][0], val['AD'][0][0],val['AD'][0][1], val['DP'][0])])
for x in sorted(myList, key = lambda x: (x[0], int(x[1]))):
    VCF_corrA.write('%s\n' %('\t'.join(str(y) for y in x)))
    
print('\tnumber of SNPs in haplotype A:', len(vcfA_dict)) #.pass.snp.vcf
print('\t\tnumber of corrected SNPs (solved + unphased (-amb1,-amb2)):',len(vcfA_corr_dict))
print('\t\tnumber of unsolved:', len(unsolvedA))

######## to obtain and print corrected VCF from haplotype B ##################################
print('\n\thaplotype B')
res= VCF_corr(vcfB_dict,vcfA_dict,randomDict,'B', coord) #obtain the corrected VCF; WARNING: the order of the dic is important!! ####new
vcfB_corr_dict=res[0]
randomDict=res[1]
unsolvedB=res[2]
unphasedB=res[3]

tag=os.path.abspath(args.VCF_B)
m=re.search('(.*)\.pass\.snp\.vcf', tag)
tag=m.group(1)
unsolvedBout=open('%s.unsolved_amb%s.bed' % (tag, ambiguity), 'w')
unphasedBout=open('%s.unphasedVariants_amb%s.bed' % (tag, ambiguity), 'w')

for x in sorted(unsolvedB, key=lambda x: (x[0], int(x[1]))): #to print SORTED chr and position of unsolved cases (see VCF_corr function)
    unsolvedBout.write('%s\t%s\t%s\n' % (x[0],int(x[1])-1,x[1]))

for x in sorted(unphasedB, key=lambda x: (x[0], int(x[1]))): #to print SORTED chr and position of unphased cases (see VCF_corr function)
    unphasedBout.write('%s\t%s\t%s\n' % (x[0],int(x[1])-1,x[1]))
    
VCF_corrB= open(VCF_corrB, 'w')
parse_vcf_header(VCF_rawB,VCF_corrB,ambiguity) #parse and write header in VCF corr

### to print the sorted vcf file:
myList=[]
for key, val in vcfB_corr_dict.items():
    if val['GT'][0] ==1:#to append only positions with GT=1
        myList.append([val['chr'][0], val['position'][0],'.',val['reference'][0], val['alternative'][0],'.','.','.','GT:AD:DP\t%s:%s,%s:%s' % (val['GT'][0], val['AD'][0][0],val['AD'][0][1], val['DP'][0])])
for x in sorted(myList, key = lambda x: (x[0], int(x[1]))):
    VCF_corrB.write('%s\n' %('\t'.join(str(y) for y in x)))

print('\tnumber of SNPs in haplotype B:' , len(vcfB_dict)) #.pass.snp.vcf
print('\t\tnumber of corrected SNPs (solved + unphased (-amb1,-amb2)):',len(vcfB_corr_dict))
print('\t\tnumber of unsolved:', len(unsolvedB))
