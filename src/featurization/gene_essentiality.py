from collections import defaultdict
import numpy as np
import os
import pandas as pd
from scipy import stats

import utils
import metadata

oe_transcript_local_path = utils.get_absolute_path_from_here("external_data/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz")
oe_gene_local_path = utils.get_absolute_path_from_here("external_data/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz")

oe_transcript_url = "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz"
oe_gene_url = "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"

gene_oe_percentiles = None
transcript_oe_percentiles = None

def init(reference_data): 
    # Download data if not already present
    utils.download_data_if_necessary(oe_transcript_local_path, oe_transcript_url)
    utils.download_data_if_necessary(oe_gene_local_path, oe_gene_url)

    # Get expected gene and transcript LOF counts
    gene_df = pd.read_csv(oe_gene_local_path, compression='gzip', sep='\t')
    transcript_df = pd.read_csv(oe_transcript_local_path, compression='gzip', sep='\t')
   
    gene_exp_counts = pd.Series(gene_df.exp_lof.values, index=gene_df.gene, dtype=float).to_dict()
    transcript_exp_counts = pd.Series(transcript_df.exp_lof.values, index=transcript_df.transcript, dtype=float).to_dict()

    # Get observed counts
    gene_obs_counts, transcript_obs_counts = get_observed_counts(reference_data)
    
    # Calculate observed / expected ratios
    gene_oe_ratios = {}
    for g in gene_exp_counts:
        if not np.isnan(gene_exp_counts[g]):
            gene_oe_ratios[g] = float(gene_obs_counts[g]) / float(gene_exp_counts[g])

    transcript_oe_ratios = {}
    for t in transcript_exp_counts:
        if not np.isnan(transcript_exp_counts[t]):
            transcript_oe_ratios[t] = float(transcript_obs_counts[t]) / float(transcript_exp_counts[t])
 
    # Calculate percentiles
    global gene_oe_percentiles
    global transcript_oe_percentiles
    gene_oe_percentiles = calculate_percentiles(gene_oe_ratios)
    transcript_oe_percentiles = calculate_percentiles(transcript_oe_ratios)

def get_observed_counts(reference_data):
    gene_counts = defaultdict(int)
    transcript_counts = defaultdict(int)

    for row in reference_data:
        if metadata.get_label(row) == 1:
            continue
        ac = utils.get_value_from_info(row[-1], "AC", int)
        
        seen_genes = set()
        for g in utils.get_gene_symbols_from_coding_changes(row):
            if g not in seen_genes:
                seen_genes.add(g)
                gene_counts[g] += ac

        seen_transcripts = set()
        for t in utils.get_transcript_ids_from_coding_changes(row):
            t = t.split('.')[0] # normalize to be consistent with gnomAD oe file
            if t not in seen_transcripts:
                seen_transcripts.add(t)
                transcript_counts[t] += ac

    return gene_counts, transcript_counts

def calculate_percentiles(oe_ratios_dict):
    keys, values = oe_ratios_dict.keys(), oe_ratios_dict.values()
    percentiles = stats.rankdata(list(values), method='average') / len(values) * 100.0
    return dict(zip(keys, percentiles))

# This value is averaged over all transcripts that the variant overlaps.
def get_average_transcript_oe_lof(row):
    oes = []
    for t in utils.get_transcript_ids_from_coding_changes(row):
        t = t.split('.')[0] # normalize to be consistent with gnomAD oe file
        oes.append(transcript_oe_percentiles[t] if t in transcript_oe_percentiles else 50.0)
    return np.mean(oes)

# This value is averaged over all genes that the variant overlaps.
def get_average_gene_oe_lof(row):
    oes = []
    for g in utils.get_gene_symbols_from_coding_changes(row):
        oes.append(gene_oe_percentiles[g] if g in gene_oe_percentiles else 50.0)
    return np.mean(oes)
