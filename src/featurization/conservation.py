import numpy as np
import os
import pyBigWig
import utils

phylop_local_path = utils.get_absolute_path_from_here("external_data/hg38.phyloP100way.bw")
phastcons_local_path = utils.get_absolute_path_from_here("external_data/hg38.phastCons100way.bw")

phylop_url = "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw"
phastcons_url = "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons100way/hg38.phastCons100way.bw"

phylop_bw = None
phastcons_bw = None
mean_phylop_cache = {}
mean_phastcons_cache = {}
max_phylop_cache = {}
max_phastcons_cache = {}
verbose = False

def init_bws():
    # Download data if it hasn't been yet
    utils.download_data_if_necessary(phylop_local_path, phylop_url)
    utils.download_data_if_necessary(phastcons_local_path, phastcons_url)
    
    # Create objects to read BigWig files
    global phylop_bw
    global phastcons_bw
    phylop_bw = pyBigWig.open(phylop_local_path)
    phastcons_bw = pyBigWig.open(phastcons_local_path)

def get_region_mean_phylop(chrom, start, end):
    global mean_phylop_cache
    k = (chrom, start, end)
    if k not in mean_phylop_cache:
        mean = phylop_bw.stats('chr' + chrom, start, end)[0]
        mean_phylop_cache[k] = mean if mean is not None else 0
    return mean_phylop_cache[k]

def get_region_mean_phastcons(chrom, start, end):
    global mean_phastcons_cache
    k = (chrom, start, end)
    if k not in mean_phastcons_cache:
        mean = phastcons_bw.stats('chr' + chrom, start, end)[0]
        mean_phastcons_cache[k] = mean if mean is not None else 0
    return mean_phastcons_cache[k]

def get_region_max_phylop(chrom, start, end):
    global max_phylop_cache
    k = (chrom, start, end)
    if k not in max_phylop_cache:
        max_ = phylop_bw.stats('chr' + chrom, start, end, type='max')[0]
        max_phylop_cache[k] = max_ if max_ is not None else 0
    return max_phylop_cache[k]

def get_region_max_phastcons(chrom, start, end):
    global max_phastcons_cache
    k = (chrom, start, end)
    if k not in max_phastcons_cache:
        max_ = phastcons_bw.stats('chr' + chrom, start, end, type='max')[0]
        max_phastcons_cache[k] = max_ if max_ is not None else 0
    return max_phastcons_cache[k]

# Averaged over all transcripts
def get_mean_phylop_on_overlapped_exon(row):
    vt_loc = int(row[4]) - 1
    mean_phylops = []
    for t in utils.get_transcripts_from_coding_changes(row):
        exon_start, exon_end = t.get_exon_overlapped(vt_loc)
        mean_phylop = get_region_mean_phylop(t.chrom, exon_start, exon_end)
        mean_phylops.append(mean_phylop)
    return np.mean(mean_phylops)

# Averaged over all transcripts
def get_mean_phastcons_on_overlapped_exon(row):
    vt_loc = int(row[4]) - 1
    all_mean_phastcons = []
    for t in utils.get_transcripts_from_coding_changes(row):
        exon_start, exon_end = t.get_exon_overlapped(vt_loc)
        mean_phastcons = get_region_mean_phastcons(t.chrom, exon_start, exon_end)
        all_mean_phastcons.append(mean_phastcons)
    return np.mean(all_mean_phastcons)

def get_mean_over_regions(coding_regions, means):
    assert(len(coding_regions) == len(means))
    total, lengths = 0, 0
    for region, mean in zip(coding_regions, means):
        total += (region[1] - region[0]) * mean
        lengths += region[1] - region[0]
    # Just in case there are no bp in coding regions, return 0. This should only happen if the
    # entire protein is truncated.
    if lengths == 0:
        if verbose:
            print("""Setting mean phylop/phastcons to 0 because coding region is empty. This should
                    only happen when the entire protein is being truncated.""")
        return 0
    return total / lengths

def get_mean_phylop(row, upstream=True):
    vt_loc = int(row[4]) - 1
    mean_phylops = []
    for t in utils.get_transcripts_from_coding_changes(row):
        regions = t.get_upstream_coding_regions(vt_loc) if upstream else t.get_downstream_coding_regions(vt_loc)
        phylop_per_region = [get_region_mean_phylop(t.chrom, start, end) for start, end in regions]
        mean_phylops.append(get_mean_over_regions(regions, phylop_per_region))
    return np.mean(mean_phylops)

def get_mean_phastcons(row, upstream=True):
    vt_loc = int(row[4]) - 1
    mean_phastcons = []
    for t in utils.get_transcripts_from_coding_changes(row):
        regions = t.get_upstream_coding_regions(vt_loc) if upstream else t.get_downstream_coding_regions(vt_loc)
        phastcons_per_region = [get_region_mean_phastcons(t.chrom, start, end) for start, end in regions]
        mean_phastcons.append(get_mean_over_regions(regions, phastcons_per_region))
    return np.mean(mean_phastcons)

def get_mean_phylop_diff(row):
    vt_loc = int(row[4]) - 1
    diffs = []
    for t in utils.get_transcripts_from_coding_changes(row):
        upstream_regions = t.get_upstream_coding_regions(vt_loc)
        downstream_regions = t.get_downstream_coding_regions(vt_loc)
        upstream_phylop = get_mean_over_regions(upstream_regions,
                [get_region_mean_phylop(t.chrom, start, end) for start, end in upstream_regions])
        downstream_phylop = get_mean_over_regions(downstream_regions,
                [get_region_mean_phylop(t.chrom, start, end) for start, end in downstream_regions])
        diffs.append(upstream_phylop - downstream_phylop)
    return np.mean(diffs)

def get_mean_phastcons_diff(row):
    vt_loc = int(row[4]) - 1
    diffs = []
    for t in utils.get_transcripts_from_coding_changes(row):
        upstream_regions = t.get_upstream_coding_regions(vt_loc)
        downstream_regions = t.get_downstream_coding_regions(vt_loc)
        upstream_phastcons = get_mean_over_regions(upstream_regions,
                [get_region_mean_phastcons(t.chrom, start, end) for start, end in upstream_regions])
        downstream_phastcons = get_mean_over_regions(downstream_regions,
                [get_region_mean_phastcons(t.chrom, start, end) for start, end in downstream_regions])
        diffs.append(upstream_phastcons - downstream_phastcons)
    return np.mean(diffs)

def get_max_phylop_diff(row):
    vt_loc = int(row[4]) - 1
    diffs = []
    for t in utils.get_transcripts_from_coding_changes(row):
        upstream_regions = t.get_upstream_coding_regions(vt_loc)
        downstream_regions = t.get_downstream_coding_regions(vt_loc)
        upstream_max, downstream_max = 0, 0
        if len(upstream_regions) > 0:
            upstream_max = max([get_region_max_phylop(t.chrom, start, end)
                    for start, end in upstream_regions])
        if len(downstream_regions) > 0:
            downstream_max = max([get_region_max_phylop(t.chrom, start, end)
                    for start, end in downstream_regions])
        diffs.append(upstream_max - downstream_max)
    return np.mean(diffs)

def get_max_phastcons_diff(row):
    vt_loc = int(row[4]) - 1
    diffs = []
    for t in utils.get_transcripts_from_coding_changes(row):
        upstream_regions = t.get_upstream_coding_regions(vt_loc)
        downstream_regions = t.get_downstream_coding_regions(vt_loc)    
        upstream_max, downstream_max = 0, 0
        if len(upstream_regions) > 0:
            upstream_max = max([get_region_max_phastcons(t.chrom, start, end) 
                    for start, end in upstream_regions])
        if len(downstream_regions) > 0:
            downstream_max = max([get_region_max_phastcons(t.chrom, start, end)
                    for start, end in downstream_regions])
        diffs.append(upstream_max - downstream_max)
    return np.mean(diffs)
