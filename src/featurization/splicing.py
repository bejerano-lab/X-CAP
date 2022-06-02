from collections import defaultdict
import utils

# global variables
gene_to_transcript_map = None
transcript_to_gene_map = None

def init(attributes_file):
    global gene_to_transcript_map
    global transcript_to_gene_map
    
    transcript_to_gene_map = {}
    gene_to_transcript_map = defaultdict(list)
    
    for line in open(attributes_file, 'r'):
        row = line.strip().split('\t')
        gene = row[0]
        transcript = row[4]
        transcript_to_gene_map[transcript] = gene
        gene_to_transcript_map[gene].append(transcript)

def can_be_spliced_out(row):
    affected_transcripts = set(utils.get_transcript_ids_from_coding_changes(row))
    affected_genes = set([transcript_to_gene_map[t] for t in affected_transcripts])
    
    potential_transcripts = set()
    for g in affected_genes:
        potential_transcripts.update(gene_to_transcript_map[g])
    
    vt_loc = int(row[4]) - 1
    for tid in potential_transcripts:
        t = utils.get_transcript(tid, row[3])
        if t is not None and not t.contains_coding_variant(vt_loc):
            return int(True)
    
    return int(False) 
