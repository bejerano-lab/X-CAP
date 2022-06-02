import utils
from utils import Variant
import metadata
from collections import defaultdict

# global variables (specific to train variants)
transcript_counts = defaultdict(lambda: defaultdict(int))
transcript_exon_counts = defaultdict(lambda: defaultdict(int))

def init(reference_data):
    for row in reference_data:
        label = bool(metadata.get_label(row))
        variant_location = int(row[4]) - 1
        for transcript in utils.get_transcripts_from_coding_changes(row):
            exon_num = transcript.get_overlapped_exon_number(variant_location)
            transcript_counts[(transcript.tid, transcript.chrom)][label] += 1
            transcript_exon_counts[(transcript.tid, transcript.chrom, exon_num)][label] += 1

def on_monoclass_pathogenic_transcript(row, is_train=False):
    is_benign = not bool(metadata.get_label(row))
    for t in utils.get_transcripts_from_coding_changes(row):
        num_benign = transcript_counts[(t.tid, t.chrom)][False]
        num_path = transcript_counts[(t.tid, t.chrom)][True]
        if is_train:
            # In training, discard impact of current variant.
            if is_benign:
                num_benign -= 1
            else:
                num_path -= 1
        if num_benign == 0 and num_path > 0:
            return 1
    return 0

def on_monoclass_pathogenic_exon(row, is_train=False):
    is_benign = not bool(metadata.get_label(row))
    variant_location = int(row[4]) - 1
    for t in utils.get_transcripts_from_coding_changes(row):
        exon_num = t.get_overlapped_exon_number(variant_location)
        num_benign = transcript_exon_counts[(t.tid, t.chrom, exon_num)][False]
        num_path = transcript_exon_counts[(t.tid, t.chrom, exon_num)][True]
        if is_train:
            # In training, discard impact of current variant.
            if is_benign:
                num_benign -= 1
            else:
                num_path -= 1
        if num_benign == 0 and num_path > 0:
            return 1
    return 0
