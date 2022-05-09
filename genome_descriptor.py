#!/usr/bin/env python3
# -*- coding=utf-8 -*-

import sys, gzip, pprint, re, argparse
import numpy as np
import pandas as pd

GEN_BIN = {
#     "./.":"?",
#     ".|.":"?",
    "./.":"0",
    ".|.":"0",
    "0/0":"0",
    "0|0":"0",
    "0/1":"1",
    "0|1":"1",
    "1/0":"1",
    "1|0":"1",
    "1/1":"2",
    "1|1":"2",
}

def site_location(site: int, start: int, end: int) -> bool:
    """snp位点是否在某个区域内，返回bool值"""
    return True if start <= site <= end else False

def gentoye2genobin(genotypes):
    return [GEN_BIN[genotype] for genotype in genotypes]

def vcf_info_parse(line):
    vcf_info_list = line.strip().split('\t')
    chorm = vcf_info_list[0]
    position = vcf_info_list[1]
    sample_genotype = gentoye2genobin([varients.split(':')[0] for varients in vcf_info_list[9:]])
    return [chorm ,position, np.array(sample_genotype,dtype='int8')]

def gff_info_parse(line):
    chrom = 'chr'+ str(line[0]).rjust(2,'0') # gff文件中的染色体标识与vcf文件中不一致，需要考虑比较妥善的处理方式
    start = line[3]
    end = line[4]
    geneid = re.match('ID=gene:([0-9A-Za-z]+);',line[8]).group(1)
    return [chrom, start, end, geneid]    


def main():
    parser = argparse.ArgumentParser(description='Descriptor for Genome.')
    parser.add_argument('--gzvcf',  help='The gzvcf file', required=True)
    parser.add_argument('--gff3',  help='The gff3 file', required=True)
    parser.add_argument('--outfile', help='The output file with csv format', required=True)
    args = parser.parse_args()
    
    gfffile = args.gff3
    vcffile = args.gzvcf
    meta_info = []
    header_info = []
    descriptor = {}
    with open(gfffile) as gffs:
        for gff_records in gffs:
            if not gff_records.startswith("#"):
                gff_record = gff_records.strip().split('\t')
                biotype = gff_record[2]
                if biotype == 'gene':
                    chrom_gff, star, end ,geneid = gff_info_parse(gff_record)
                    descriptor[geneid] = np.array([0])
                    with gzip.open(vcffile, mode='rt') as vcf:
                        for line in vcf:
                            if not line.startswith("#"):
                                chrom_vcf, position, genotypes = vcf_info_parse(line)
                                if chrom_vcf == chrom_gff: # 如果变异位置是某个基因中则记录，超过这个基因的结尾则退出当前循环
                                    if position < star:  
                                        continue
                                    elif site_location(position,star,end):
                                        descriptor[geneid] = descriptor[geneid] + genotypes #利用array的广播机制记录vcf的杂合性
                                    else:
                                        break
                            elif line.startswith("##"):
                                meta_info.append(line.strip())
                            else:
                                header_info.append(line.strip().split('\t'))
    pd.DataFrame(descriptor,index=header_info[9:]).T.to_csv(args.outfile)

if __name__ == "__main__":
    main()

with gzip.open(vcffile, mode='rt') as vcf:
    for line in vcf:
        if not line.startswith("#"):
            chrom_vcf, position, genotypes = vcf_info_parse(line)
            if chrom_vcf == chrom_gff: # 如果变异位置是某个基因中则记录，超过这个基因的结尾则退出当前循环
                if position < star:  
                    continue
                elif site_location(position,star,end):
                    descriptor[geneid] = descriptor[geneid] + genotypes #利用array的广播机制记录vcf的杂合性
                else:
                    break
        elif line.startswith("##"):
            meta_info.append(line.strip())
        else:
            header_info.append(line.strip().split('\t'))