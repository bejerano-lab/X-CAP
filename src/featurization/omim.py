import utils

omim_filepath = utils.get_absolute_path_from_here("external_data/genemap2.txt")

# global variables
recessive_genes = None
dominant_genes = None

def init():
    global recessive_genes
    global dominant_genes
    
    header_map = None
    recessive_genes = set()
    dominant_genes = set()

    for line in open(omim_filepath):
        row = line.strip().split('\t')
        if line.startswith("#"):
            if "Gene Symbols" in line:
                header_map = {hdr: i for i, hdr in enumerate(row)}

        elif max(header_map["Gene Symbols"], header_map["Phenotypes"]) < len(row):
            gene_symbols = row[header_map["Gene Symbols"]]
            phenotypes = row[header_map["Phenotypes"]]
            recessive = "recessive" in phenotypes
            dominant = "dominant" in phenotypes
            for gene in gene_symbols.split(','):
                gene = gene.strip()
                if recessive:
                    recessive_genes.add(gene)
                if dominant:
                    dominant_genes.add(gene)
 
def on_recessive_gene(row):
    on_recessive = any([g in recessive_genes for g in utils.get_gene_symbols_from_coding_changes(row)])
    return int(on_recessive)

def on_dominant_gene(row):
    on_dominant = any([g in dominant_genes for g in utils.get_gene_symbols_from_coding_changes(row)])
    return int(on_dominant)
