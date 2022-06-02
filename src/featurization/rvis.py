import numpy as np
import utils

rvis_path = utils.get_absolute_path_from_here("external_data/RVIS_Unpublished_ExACv2_March2017.txt")

# global variables
gene_to_rvis = None

def init():
    global gene_to_rvis
    
    rvis_file = open(rvis_path)
    
    header = next(rvis_file).strip().split('\t')
    gene_index = header.index("CCDSr20")
    rvis_index = header.index("%RVIS[pop_maf_0.05%(any)]")
    
    gene_to_rvis = {}
    for line in rvis_file:
        row = line.strip().split('\t')
        gene = row[gene_index]
        try:
            gene_to_rvis[gene] = float(row[rvis_index])
        except ValueError:
            continue
    
def get_gene_rvis(gene):
    if gene in gene_to_rvis:
        return gene_to_rvis[gene]
    else:
        return 50.0

# This value is averaged over all genes that the variant overlaps.
# Lower values mean more intolerant.
def get_average_rvis(row):
    percentiles = []
    for gene in utils.get_gene_symbols_from_coding_changes(row):
        percentile = get_gene_rvis(gene)
        percentiles.append(percentile) 
    return np.mean(percentiles)
