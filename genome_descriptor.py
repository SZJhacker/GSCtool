#!/usr/bin/env python3
# -*- coding=utf-8 -*-

"""
The genic descriptor characterizes the genome in a time-series manner. The script counts the number of
single nucleotide polymorphisms(SNPs) in per geneid of VCF files depend on
GFF3 record.
"""

__author__      = "Zijie Shen"
__version__     = "1.0"
__email__       = "shenzijie2013@163.com"
__date__        = "2022-05-09"

import gzip, re, argparse
import numpy as np
import pandas as pd

GEN_BIN = {
    "./.":0,
    ".|.":0,
    "0/0":0,
    "0|0":0,
    "0/1":1,
    "0|1":1,
    "1/0":1,
    "1|0":1,
    "1/1":2,
    "1|1":2,
}

def site_location(site: int, start: int, end: int) -> bool:
    """snp位点是否在某个区域内，返回bool值
    Determine whether the snp site is in a certain region, return bool values"""
    return True if start <= site <= end else False

def gentoye2genobin(genotypes):
    """converted the snps into genotype"""
    return [GEN_BIN.get(genotype,0) for genotype in genotypes]


def vcf_info_parse(line):
    vcf_info_list = line.strip().split('\t')
    chorm = vcf_info_list[0]
    position = vcf_info_list[1]
    sample_genotype = gentoye2genobin([varients.split(':')[0] for varients in vcf_info_list[9:]])
    return [chorm,position, np.array(sample_genotype,dtype=np.int8)]

def gff_info_parse(line):
    chrom = line[0]
    start = line[3]
    end = line[4]
    geneid_pre = re.search('ID=([0-9A-Za-z:]+);',line[8]).group(1)
    geneid = re.sub('gene:','',geneid_pre,re.I)
    return [chrom, start, end, geneid]    


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gzvcf',  help='The gzvcf file of individual or population', required=True)
    parser.add_argument('--gff3',  help='The gff3 file which corresponds to the spieces of the gzvcf', required=True)
    parser.add_argument('--outfile', help='The name of output file with csv format, which records the features of genome', required=True)
    args = parser.parse_args()
    
    gfffile = args.gff3
    vcffile = args.gzvcf


    meta_info = []
    header_info = []
    descriptor = {}
    anchor = False
    with open(gfffile) as gffs, gzip.open(vcffile, mode='rt') as vcf:
        for gff_records in gffs:
            gff_record = gff_records.strip().split('\t')
            biotype = gff_record[2:3] # 获得索引为2的值，越界不报错(注释行可能没有9个元素)
            if biotype and biotype[0] == 'gene': 
                chrom_gff, star, end ,geneid = gff_info_parse(gff_record)
                descriptor[geneid] = np.array([0])
                if anchor and site_location(position,star,end): 
                    descriptor[geneid] = descriptor[geneid] + genotypes #利用array的广播机制记录vcf的杂合性
                    anchor = False
                for line in vcf:
                    if not line.startswith("#"):
                        chrom_vcf, position, genotypes = vcf_info_parse(line)
                        if chrom_vcf == chrom_gff:
                            if position < star:  
                                continue
                            elif site_location(position,star,end):
                                descriptor[geneid] = descriptor[geneid] + genotypes #利用array的广播机制记录vcf的杂合性
                            else:
                                anchor=True
                                break
                        else:
                            anchor=True # 染色体不一致的时候，vcf迭代到下一条染色体，跳出vcf文件的循环后应该进行一次vcf和gff文件的操作
                            break
                    elif line.startswith("##"):
                        meta_info.append(line.strip())
                    else:
                        header_info.append(line.strip())

    pd.DataFrame(descriptor,index=header_info[0].split('\t')[9:]).to_csv(args.outfile)

if __name__ == "__main__":
    main()
