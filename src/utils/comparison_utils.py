import itertools
import os
from annovar_utils import *
from variant_utils import Variant

xcap_train_path = "/cluster/u/rrastogi/ECNN/results/d_original/xcap/train_set/filtered_train.tsv"
_1KG_dir = "/cluster/u/rrastogi/data/1KG/phase1/hg38"
ExAC_file = "/cluster/u/rrastogi/data/ExAC_0.3/ExAC.0.3.GRCH38.condensed.vcf"
hgmd_2016_path = "/cluster/data/labResources/medgen/hgmd/HGMD_PRO_2016.2/HGMD_PRO_2016.2_hg38.vcf"

def get_xcap_train_variants():
    positive_variants, negative_variants = set(), set()
    for line in open(xcap_train_path):
        v = get_variant_from_line(line)
        positive = get_value_from_line(line, "LABEL", lambda x: bool(int(x)))
        if positive:
            positive_variants.add(v)
        else:
            negative_variants.add(v)
    return positive_variants, negative_variants

def get_competitors_positive_variants():
    positive_variants = set()
    for line in open(hgmd_2016_path):
        if line.startswith("#"):
            continue
        row = line.strip().split('\t')
        v = Variant(row[0], row[1], row[3], row[4])
        positive_variants.add(v)
    return positive_variants

# This returns a set of SNPs within phase1 1KG and ExAC -- not just stopgains
def get_competitors_negative_variants():
    negative_variants = set()
    chromosomes = list(map(str, range(1, 23))) + ['X'] # no chrom Y data in phase1 1KG
    for chrom in chromosomes:
        filepath = os.path.join(_1KG_dir, "chr{}.hg38.vcf".format(chrom))
        for line in open(filepath):
            if line.startswith("#"):
                continue
            row = line.strip().split('\t')
            v = Variant(row[0], row[1], row[3], row[4])
            if len(v.ref) == 1 and len(v.alt) == 1:
                negative_variants.add(v)

    for line in open(ExAC_file):
        row = line.strip().split('\t')
        refs, alts = row[2].split(','), row[3].split(',')
        for ref, alt in itertools.product(refs, alts):
            v = Variant(row[0], row[1], ref, alt)
            if len(v.ref) == 1 and len(v.alt) == 1:
                negative_variants.add(v)

    return negative_variants

def get_all_potential_train_variants():
    print("getting X-CAP variants")
    xcap_train_positive_variants, xcap_train_negative_variants = get_xcap_train_variants()
    print("getting others' positive variants")
    comp_positive_variants = get_competitors_positive_variants()
    print("getting others' negative variants")
    comp_negative_variants = get_competitors_negative_variants()
    return xcap_train_positive_variants | xcap_train_negative_variants | comp_positive_variants | comp_negative_variants
