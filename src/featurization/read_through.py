import itertools
import utils

complement_map = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C"
}

stop_codons = ["TAG", "TAA", "TGA"]

kmers_cache = {}

def generate_kmers(num_context):
    if num_context not in kmers_cache:
        context = [''.join(p) for p in itertools.product('ACTG', repeat=num_context)]
        kmers = [f'{a}{b}' for a in stop_codons for b in context]
        kmers_cache[num_context] = kmers
    return kmers_cache[num_context]

def get_stop_codon(row, num_context=0):
    kmers = generate_kmers(num_context)
    one_hot_encoding = [0 for _ in range(len(kmers))]

    variant_loc = int(row[4]) - 1
    for t in utils.get_transcripts_from_coding_changes(row):
        variant_index = t.determine_mrna_index(variant_loc)
        modified_index = (variant_index - t.cds_start_index) % 3
        codon_start = variant_index - modified_index
        codon_end = codon_start + 3
        
        alt = row[7] if t.strand == "+" else complement_map[row[7]]
        seq = t.mrna[codon_start: codon_start + modified_index] + alt + \
              t.mrna[codon_start + modified_index + 1: codon_end] + \
              t.mrna[codon_end: codon_end + num_context]

        if len(seq) == 3 + num_context:
            one_hot_encoding[kmers.index(seq)] = 1
            break # once we have one, break early

    return kmers, one_hot_encoding
