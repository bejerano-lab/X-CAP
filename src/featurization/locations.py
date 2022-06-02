import numpy as np
import utils
from utils import Transcript

import sys
sys.path.append(utils.get_absolute_path_from_here("../utils"))
from variant_utils import Variant

def get_distance_from_cds_start(row):
    vt_loc = int(row[4]) - 1
    distances = []
    for t in utils.get_transcripts_from_coding_changes(row):
        distance = t.get_location_within_transcript(vt_loc)
        distances.append(distance)
    return np.mean(distances)

def get_distance_from_cds_end(row):
    vt_loc = int(row[4]) - 1
    distances = []
    for t in utils.get_transcripts_from_coding_changes(row):
        distance = t.get_coding_length() - t.get_location_within_transcript(vt_loc)
        distances.append(distance)
    return np.mean(distances)

def calculate_relative_location_in_cds(t, vt_loc):
    return t.get_location_within_transcript(vt_loc) / t.get_coding_length()

def get_relative_location_in_cds(row):
    vt_loc = int(row[4]) - 1 
    relative_locations = []
    for t in utils.get_transcripts_from_coding_changes(row):
        relative_locations.append(calculate_relative_location_in_cds(t, vt_loc)) 
    return np.mean(relative_locations)

def get_overlapped_exon_num(row):
    vt_loc = int(row[4]) - 1
    overlapped_exon_nums = []
    for t in utils.get_transcripts_from_coding_changes(row):
        overlapped_exon_nums.append(t.get_overlapped_exon_number(vt_loc))
    return np.mean(overlapped_exon_nums)

def get_overlapped_exon_length(row):
    vt_loc = int(row[4]) - 1
    overlapped_exon_lengths = []
    for t in utils.get_transcripts_from_coding_changes(row):
        overlapped_exon_lengths.append(t.get_overlapped_exon_length(vt_loc))
    return np.mean(overlapped_exon_lengths)

def get_distance_from_exon_start(row):
    vt_loc = int(row[4]) - 1
    distances = []
    for t in utils.get_transcripts_from_coding_changes(row):
        distance = t.get_location_within_exon(vt_loc)
        distances.append(t.get_location_within_exon(vt_loc))
    return np.mean(distances)

def get_distance_from_exon_end(row):
    vt_loc = int(row[4]) - 1
    distances = []
    for t in utils.get_transcripts_from_coding_changes(row):
        distance = t.get_overlapped_exon_length(vt_loc) - t.get_location_within_exon(vt_loc) - 1
        distances.append(distance)
    return np.mean(distances)

def calculate_relative_location_in_exon(t, vt_loc):
    return t.get_location_within_exon(vt_loc) / t.get_overlapped_exon_length(vt_loc)

def get_relative_location_in_exon(row):
    vt_loc = int(row[4]) - 1
    relative_locations = []
    for t in utils.get_transcripts_from_coding_changes(row):
        relative_locations.append(calculate_relative_location_in_exon(t, vt_loc))
    return np.mean(relative_locations)

def get_percentage_of_transcripts_with_nmd(row):
    vt_loc = int(row[4]) - 1 
    nmd_is_triggered = []
    for t in utils.get_transcripts_from_coding_changes(row):
        nmd_is_triggered.append(int(t.triggers_nmd(vt_loc)))
    return np.mean(nmd_is_triggered)

def get_distance_from_last_exon_exon_junction(row):
    vt_loc = int(row[4]) - 1 
    return np.mean([t.get_distance_from_last_exon_exon_junction(vt_loc)
                    for t in utils.get_transcripts_from_coding_changes(row)])

def on_autosomal_chromosome(row):
    v = Variant(row[3], row[4], row[6], row[7])
    return int(v.chrom not in ["X", "Y"])

def on_X_chromosome(row):
    v = Variant(row[3], row[4], row[6], row[7])
    return int(v.chrom == "X")

def on_Y_chromosome(row):
    v = Variant(row[3], row[4], row[6], row[7])
    return int(v.chrom == "Y")
