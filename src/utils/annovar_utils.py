from variant_utils import Variant

def get_value_from_line(line, key, dtype=str):
    row = line.strip().split('\t')
    return get_value_from_info(row[-1], key, dtype)

def get_value_from_info(info, key, dtype=str):
    for it in info.split(';'):
        if '=' not in it:
            continue
        k, v = it.split('=')
        if k == key:
            return dtype(v)

def get_variant_from_line(line):
    row = line.strip().split('\t')
    return Variant(row[3], row[4], row[6], row[7])

def get_transcripts_from_coding_changes(changes):
    transcripts = []
    for change in changes.split(','):
        change = change.strip()
        if len(change) != 0:
            t = change.split(':')[1]
            transcripts.append(t)
    return transcripts

def get_gene_symbols_from_coding_changes(changes):
    genes = []
    for change in changes.split(','):
        if len(change) != 0:
            g = change.split(':')[0]
            genes.append(g)
    return genes
