import numpy as np
import utils

def get_cds_length(row):
    cds_lengths = []
    for t in utils.get_transcripts_from_coding_changes(row):
        cds_lengths.append(t.get_coding_length())
    return np.mean(cds_lengths)

def get_num_exons(row):
    num_exons = []
    for t in utils.get_transcripts_from_coding_changes(row):
        num_exons.append(t.get_num_exons())
    return np.mean(num_exons)

if __name__ == '__main__':
    main()
