import os
import subprocess

import sys
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "utils"))
from variant_utils import Variant

stop_codons = ["TAG", "TAA", "TGA"]
nucleotides = ["A", "T", "C", "G"]
complement = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C"
}

class Transcript():
    def __init__(self, tid, chrom, strand):
        self.tid = tid
        self.chrom = chrom.replace('chr', '')
        self.strand = strand
        
        self.cds_start = None
        self.cds_end = None
        self.exon_starts = None
        self.exon_ends = None
        self.coding_exon_starts = None
        self.coding_exon_ends = None
        self.mrna = None # contains UTRs as well
        self.cds_start_index = None # index of cds_start in self.mrna
        self.cds_end_index = None # index of cds_end in self.mrna (equal to 1 after stop codon)
    
    def set_cds_span(self, cds_start, cds_end):
        self.cds_start = int(cds_start)
        self.cds_end = int(cds_end)

    def set_exons(self, exon_starts, exon_ends):        
        self.exon_starts = list(map(int, split_and_discard_empty(exon_starts)))
        self.exon_ends = list(map(int, split_and_discard_empty(exon_ends)))

    # Pre-conition: set_cds_span should have already been called
    def set_sequence(self, sequence):
        self.mrna = sequence

        real_cds_start = self.cds_start if self.strand == "+" else self.cds_end - 1
        self.cds_start_index = self.determine_mrna_index(real_cds_start)

        real_cds_end = self.cds_end - 1 if self.strand == "+" else self.cds_start
        self.cds_end_index = self.determine_mrna_index(real_cds_end) + 1

    def determine_mrna_index(self, location):
        mrna_index = 0
        for start, end in zip(self.exon_starts, self.exon_ends):
            if self.overlaps_region(location, start, end):
                mrna_index += (location - start) if self.strand == "+" else (end - 1 - location)
            if self.downstream_region(location, start, end):
                mrna_index += end - start
        return mrna_index

    def get_location_from_mrna_index(self, mrna_index):
        remaining_len = mrna_index

        exon_starts = list(self.exon_starts) if self.strand == "+" else list(reversed(self.exon_starts))
        exon_ends = list(self.exon_ends) if self.strand == "+" else list(reversed(self.exon_ends))

        for start, end in zip(exon_starts, exon_ends):
            if remaining_len >= end - start:
                remaining_len -= end - start
                continue
            else:
                location = start + remaining_len if self.strand == "+" else end - 1 - remaining_len
                return location

    def get_mrna(self):
        return self.mrna

    def get_potential_stopgains(self):
        variants = set()

        for i in range(self.cds_start_index, self.cds_end_index, 3):
            codon = self.mrna[i: i + 3]
            if i + 3 >= self.cds_end_index or codon in stop_codons:
                break

            for offset in range(len(codon)):
                for n in nucleotides:
                    new_codon = codon[: offset] + n + codon[offset + 1: ]
                    if new_codon in stop_codons:
                        location = self.get_location_from_mrna_index(i + offset) + 1 # +1 to convert to 1-index
                        ref = codon[offset] if self.strand == "+" else complement[codon[offset]]
                        alt = n if self.strand == "+" else complement[n]
                        v = Variant(self.chrom, location, ref, alt)
                        variants.add(v)

        return variants

    def get_num_exons(self):
        assert (len(self.exon_starts) == len(self.exon_ends))
        return len(self.exon_starts)
    
    # Returns a list of elements of the form (exon start, exon end). These exons are in order that
    # the transcript is transcribed. However, start < end even for transcripts on the - strand.
    def get_coding_exons(self):
        if self.coding_exon_starts and self.coding_exon_ends:
            return self.coding_exon_starts, self.coding_exon_ends

        coding_exon_starts, coding_exon_ends = [], []
        
        for start, end in zip(self.exon_starts, self.exon_ends):
            if self.cds_end < start or self.cds_start >= end:
                continue

            coding_exon_start = max(self.cds_start, start)
            coding_exon_start = min(coding_exon_start, end)
            coding_exon_starts.append(coding_exon_start)

            coding_exon_end = min(self.cds_end, end)
            coding_exon_end = max(coding_exon_end, start)
            coding_exon_ends.append(coding_exon_end)

        if self.strand == "-":
            coding_exon_starts = list(reversed(coding_exon_starts))
            coding_exon_ends = list(reversed(coding_exon_ends))

        return coding_exon_starts, coding_exon_ends
    
    # variant_location should be 0-indexed
    def contains_coding_variant(self, variant_location):
        exon_starts, exon_ends = self.get_coding_exons()
        return any([self.overlaps_region(variant_location, start, end) for start, end
                    in zip(exon_starts, exon_ends)])
    
    def get_coding_length(self):
        coding_length = 0
        exon_starts, exon_ends = self.get_coding_exons()
        for start, end in zip(exon_starts, exon_ends):
            coding_length += end - start # don't have to add 1 because end is exclusive
        return coding_length

    def overlaps_region(self, variant_location, start, end):
        return start <= variant_location < end

    def downstream_region(self, variant_location, start, end):
        if self.strand == "+":
            return variant_location >= end
        else:
            return variant_location < start

    def upstream_region(self, variant_location, start, end):
        if self.strand == "+":
            return variant_location < start
        else:
            return variant_location >= end
    
    # variant_location should be 0-indexed
    def get_location_within_transcript(self, variant_location):
        exon_starts, exon_ends = self.get_coding_exons()

        length_before = 0
        for start, end in zip(exon_starts, exon_ends):
            if self.overlaps_region(variant_location, start, end):
                length_before += (variant_location - start) if self.strand == "+" else end - 1 - variant_location
            elif self.downstream_region(variant_location, start, end):
                length_before += end - start
            else:
                break      
        return length_before

    # Assuming that a variant overlaps a transcript, returns the exon number in which the variant
    # is located.
    def get_overlapped_exon_number(self, variant_location):
        exon_starts, exon_ends = self.get_coding_exons()
        for i, (start, end) in enumerate(zip(exon_starts, exon_ends)):
            if self.overlaps_region(variant_location, start, end):
                return i + 1
        return None

    # Assuming that a variant overlaps a transcript, returns the exon_start and exon_end
    # of the coding region of the exon. The variant_location should be 0-indexed.
    def get_exon_overlapped(self, variant_location):
        exon_starts, exon_ends = self.get_coding_exons()
        for start, end in zip(exon_starts, exon_ends):
            if self.overlaps_region(variant_location, start, end):
                return start, end
        return None

    # Assuming that a variant overlaps a transcript, returns the length of the exon that is overlapped.
    # The variant location should be 0-indexed.
    def get_overlapped_exon_length(self, variant_location):
        start, end = self.get_exon_overlapped(variant_location)
        return end - start

    # variant_location should be 0-indexed
    def get_location_within_exon(self, variant_location):
        start, end = self.get_exon_overlapped(variant_location)
        return (variant_location - start) if self.strand == "+" else end - 1 - variant_location

    def remove_empty_regions(self, regions):
        ret = []
        for start, end in regions:
            if start != end:
                ret.append((start, end))
        return ret

    def get_upstream_coding_regions(self, variant_location):
        assert(self.contains_coding_variant(variant_location))

        exon_starts, exon_ends = self.get_coding_exons()
        upstream_coding_regions = [] # end of region is exclusive

        for start, end in zip(exon_starts, exon_ends):
            if self.overlaps_region(variant_location, start, end):
                region = (start, variant_location) if self.strand == "+" else (variant_location + 1, end)
                upstream_coding_regions.append(region)
            elif self.downstream_region(variant_location, start, end):
                upstream_coding_regions.append((start, end))
            else:
                break

        # remove empty regions in case variant is at very beginning of an exon
        return self.remove_empty_regions(upstream_coding_regions)

    def get_downstream_coding_regions(self, variant_location):
        assert(self.contains_coding_variant(variant_location))

        exon_starts, exon_ends = self.get_coding_exons()
        downstream_coding_regions = [] # end of region is exclusive

        for start, end in zip(exon_starts, exon_ends):
            if self.overlaps_region(variant_location, start, end):
                region = (variant_location + 1, end) if self.strand == "+" else (start, variant_location)
                downstream_coding_regions.append(region)
            if self.upstream_region(variant_location, start, end):
                downstream_coding_regions.append((start, end))

        # remove empty regions in case variant is at very beginning of an exon
        return self.remove_empty_regions(downstream_coding_regions)

    def get_last_exon_exon_junction_location(self):
        exon_starts, exon_ends = self.get_coding_exons()
        if len(exon_starts) == 1:
            return None
        return exon_ends[-2] - 1 if self.strand == "+" else exon_starts[-2]

    # variant location should be 0-indexed
    def triggers_nmd(self, variant_location):
        junction_location = self.get_last_exon_exon_junction_location()
        if junction_location is None:
            return False

        junction_index = self.get_location_within_transcript(junction_location)
        variant_index = self.get_location_within_transcript(variant_location)
        return (junction_index - variant_index) > 50

    # variant location should be 0-indexed
    def get_distance_from_last_exon_exon_junction(self, variant_location):
        junction_location = self.get_last_exon_exon_junction_location()
        if junction_location is None:
            return 0 - self.get_location_within_transcript(variant_location)

        junction_index = self.get_location_within_transcript(junction_location)
        variant_index = self.get_location_within_transcript(variant_location)
        return junction_index - variant_index

def get_transcript_ids_from_coding_changes(row):
    transcripts = []
    for change in row[2].split(','):
        if change.strip() == "":
            continue
        transcript = change.split(':')[1]
        transcripts.append(transcript)
    return transcripts

def get_gene_symbols_from_coding_changes(row):
    genes = []
    for change in row[2].split(','):
        if change.strip() == "":
            continue
        gene = change.split(':')[0]
        genes.append(gene)
    return genes

def split_and_discard_empty(arr, delimiter=','):
    splits = []
    for val in arr.split(delimiter):
        if len(val.strip()) != 0:
            splits.append(val)
    return splits

def get_value_from_info(info, key, dtype):
    for it in info.split(';'):
        if '=' not in it:
            continue
        k, v = it.split('=')
        if k == key:
            return dtype(v)

transcript_map = None

def init_transcript_map(transcripts_path, mrna_path):
    global transcript_map
    transcript_map = {}
    for line in open(transcripts_path):
        row = line.strip().split('\t')
        t = Transcript(row[1], row[2], row[3])
        t.set_cds_span(row[6], row[7])
        t.set_exons(row[9], row[10])
        transcript_map[(t.tid, t.chrom)] = t

    ensembl_id, chrom = None, None
    for line in open(mrna_path):
        if line.startswith('>'):
            ensembl_id = line.split(' ')[0].replace('>', '')
            chrom = line[line.index('chr') + 3: ].split(':')[0]
        elif (ensembl_id, chrom) in transcript_map:
            transcript_map[(ensembl_id, chrom)].set_sequence(line.strip())

# Pre-condition: init_transcript_map has been called
def get_transcript_map():
    return transcript_map

def get_transcripts_from_coding_changes(row):
    transcripts = []
    chrom = row[3].replace('chr', '')
    for change in row[2].split(','):
        if change.strip() == "":
            continue
        tid = change.split(':')[1]
        transcripts.append(transcript_map[(tid, chrom)])
    return transcripts

def get_transcript(tid, chrom):
    chrom = chrom.replace('chr', '')
    return transcript_map[(tid, chrom)] if (tid, chrom) in transcript_map else None

def get_absolute_path_from_here(local_path):
    return os.path.join(os.path.dirname(__file__), local_path)

def download_data_if_necessary(absolute_path, url): # local_path is of form "external_data/...")
    if not os.path.exists(absolute_path):
        bashCmd = "wget {}".format(url)
        external_data_dir = os.path.join(os.path.dirname(__file__), 'external_data')
        child = subprocess.Popen(bashCmd.split(), cwd=external_data_dir)
        child.wait()
        assert os.path.exists(absolute_path)
