import utils
from utils import Variant

# global variables (specific to train variants)
positive_variants = None
negative_variants = None

def init_variant_data(reference_data):
    global positive_variants
    global negative_variants
    positive_variants = set()
    negative_variants = set()
    
    for row in reference_data:
        v = Variant(row[3], row[4], row[6], row[7])
        label = get_label(row)
        if label == 1:
            positive_variants.add(v)
        else:
            negative_variants.add(v)

PAR1_X_range = (10_001, 2_781_479)
PAR1_Y_range = (10_001, 2_781_479)
PAR2_X_range = (155_701_383, 156_030_895)
PAR2_Y_range = (56_887_903, 57_217_415)

def in_range(pos, range_tuple):
    start, end = range_tuple
    return start <= pos <= end

def in_pseudoautosomal_region(v):
    if v.chrom == "X" and (in_range(v.pos, PAR1_X_range) or in_range(v.pos, PAR2_X_range)):
        return True
    if v.chrom == "Y" and (in_range(v.pos, PAR1_Y_range) or in_range(v.pos, PAR2_Y_range)):
        return True
    return False

def get_label(row):
    if "LABEL=" not in row[-1]:
        return None
    else:
        return utils.get_value_from_info(row[-1], "LABEL", int)

def get_zygosity(row):
    v = Variant(row[3], row[4], row[6], row[7])
    label = get_label(row)

    # Variants on sex chromosomes are always treated as homozygous
    if v.chrom in ["X", "Y"] and not in_pseudoautosomal_region(v):
        return 1 # homozygous
    
    # If we are given exact genotype information in the form of GT="..", use that
    if "GT=" in row[-1]:
        gt = utils.get_value_from_info(row[-1], "GT", str)
        assert(gt in ["1/1", "0/1", "1/0"])
        return int(gt == "1/1")
    
    # We know little about positive variants. Let's just see if they're in the negative set.
    # If they are, we predict that they are homozygous.
    label = get_label(row)
    assert label is not None
    if label == 1:
        return int(v in negative_variants)
    # For negative variants, we already know the truth.
    if label == 0:
        nhomalt = utils.get_value_from_info(row[-1], "NHOMALT", int)
        return int(nhomalt > 0)
