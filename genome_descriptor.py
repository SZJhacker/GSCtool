#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GSCtool (Parallel Version)

This tool characterizes genomes for machine learning applications by counting 
the number of single nucleotide polymorphisms (SNPs) within each gene region.
It uses a GFF3 file for gene annotations and a VCF file for variant calls.

This version leverages pysam for efficient VCF indexing and multiprocessing
for parallel computation, offering significant speed improvements over the
original version.
"""

__author__ = "Zijie Shen"
__version__ = "1.1"
__email__ = "shenzijie2013@163.com"
__date__ = "2025-08-07"

import pysam
import pandas as pd
import numpy as np
import re
import argparse
import multiprocessing as mp
from functools import partial


def genotype_to_score(gt):
    if gt in ["./.", ".|.", ".", ""]:
        return 0
    alleles = gt.replace('|', '/').split('/')
    if len(alleles) != 2 or not all(a.isdigit() or a == '.' for a in alleles):
        return 0
    nums = [int(a) if a != '.' else 0 for a in alleles]
    return 0 if nums == [0, 0] else (2 if nums[0] == nums[1] else 1)


def parse_gff_gene(line):
    parts = line.strip().split('\t')
    if len(parts) < 9 or parts[2] != 'gene':
        return None
    chrom, start, end, attrs = parts[0], int(parts[3]), int(parts[4]), parts[8]
    attr_dict = dict(item.split('=') for item in attrs.split(';') if '=' in item)
    gene_id = next((attr_dict.get(k, '').replace('gene:', '')
                    for k in ['ID', 'gene_id', 'Name'] if k in attr_dict), None)
    return (chrom, start, end, gene_id) if gene_id else None


def extract_sample_genotypes(record, samples):
    return [f"{gt[0]}/{gt[1]}" if gt[0] is not None else "./."
            for gt in [record.samples[s]['GT'] for s in samples]]


def load_all_genes(gff_file):
    genes = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            gene_info = parse_gff_gene(line)
            if gene_info:
                genes.append(gene_info)
    return genes


def split_genes_into_batches(genes, batch_size=1000):
    return [genes[i:i + batch_size] for i in range(0, len(genes), batch_size)]


def process_single_gene(gene_info, vcf_reader, samples):
    chrom, start, end, gene_id = gene_info
    gene_scores = np.zeros(len(samples), dtype=int)
    for record in vcf_reader.fetch(chrom, start - 1, end):
        genotypes = extract_sample_genotypes(record, samples)
        scores = [genotype_to_score(gt) for gt in genotypes]
        gene_scores += scores
    return gene_id, gene_scores


def process_gene_batch(gene_batch, vcf_file, samples):
    vcf_reader = pysam.VariantFile(vcf_file)
    results = {}
    for gene_info in gene_batch:
        try:
            gene_id, gene_scores = process_single_gene(gene_info, vcf_reader, samples)
            results[gene_id] = gene_scores
        except Exception as e:
            chrom, start, end, gene_id = gene_info
            print(f"Error processing gene {gene_id}: {e}")
            results[gene_id] = np.zeros(len(samples), dtype=int)
    vcf_reader.close()
    return results


def merge_batch_results(batch_results):
    all_results = {}
    for batch in batch_results:
        all_results.update(batch)
    return all_results


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gzvcf', required=True, help='Gzipped VCF file')
    parser.add_argument('--gff3', required=True, help='GFF3 annotation file')
    parser.add_argument('--outfile', required=True, help='Output CSV file')
    parser.add_argument('--batch_size', type=int, default=1000, help='Genes per batch (default: 1000)')
    parser.add_argument('--num_processes', type=int, default=4, help='Number of parallel processes (default: 4)')
    args = parser.parse_args()

    print("Loading samples from VCF...")
    vcf_reader = pysam.VariantFile(args.gzvcf)
    samples = list(vcf_reader.header.samples)
    vcf_reader.close()

    print("Loading gene annotations...")
    all_genes = load_all_genes(args.gff3)
    print(f"Total genes found: {len(all_genes)}")

    gene_batches = split_genes_into_batches(all_genes, args.batch_size)
    print(f"Split into {len(gene_batches)} batches with batch size {args.batch_size}")

    process_func = partial(process_gene_batch, vcf_file=args.gzvcf, samples=samples)
    print(f"Processing with {args.num_processes} processes...")
    with mp.Pool(processes=args.num_processes) as pool:
        batch_results = pool.map(process_func, gene_batches)

    print("Merging results...")
    all_results = merge_batch_results(batch_results)

    print(f"Writing to {args.outfile}")
    df = pd.DataFrame(all_results, index=samples)
    df.to_csv(args.outfile, index_label='SampleID')

    print("Done.")


if __name__ == "__main__":
    main()
