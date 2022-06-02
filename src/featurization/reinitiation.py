import numpy as np
import utils

def get_dist_to_next_start_codon(row):
    distances = []
    variant_loc = int(row[4]) - 1

    for t in utils.get_transcripts_from_coding_changes(row):
        variant_index = t.determine_mrna_index(variant_loc)
        next_codon_index = variant_index + 3 - ((variant_index - t.cds_start_index) % 3)
        try:
            distance = t.mrna[next_codon_index: t.cds_end_index].find("ATG")
        except:
            distance = t.cds_end_index - next_codon_index
        distances.append(distance)

    return np.mean(distances)
