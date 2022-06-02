import os
import pandas as pd
import subprocess
import sys

def make_directory_if_necessary(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def is_already_prepared(gencode_dir):
    fnames = set([fname for fname in os.listdir(gencode_dir)])
    return "hg38_ensGene.txt" in fnames and "hg38_ensGeneMrna.fa" in fnames

def download_gencode(annovar_dir, seqdir, gencode_tempdir):
    # 1. Download whole-genome FASTA file
    bashCmd = "perl {}/annotate_variation.pl -downdb -build hg38 seq {}".format(annovar_dir, seqdir)
    subprocess.run(bashCmd.split())

    # 2. Download Gencode transcript annotation file
    bashCmd = "perl {}/annotate_variation.pl -downdb wgEncodeGencodeBasicV33 {} -build hg38".format(annovar_dir, gencode_tempdir)
    subprocess.run(bashCmd.split())

    # 3. Download Gencode transcript attributes file
    attributes_link = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeAttrsV33.txt.gz"
    bashCmd = "wget -P {} {}".format(gencode_tempdir, attributes_link)
    subprocess.run(bashCmd.split())
    
    bashCmd = "gunzip {}/wgEncodeGencodeAttrsV33.txt.gz".format(gencode_tempdir)
    subprocess.run(bashCmd.split())

    return seqdir

def filter_transcript_set(annovar_dir, seqdir, gencode_dir, gencode_tempdir):
    transcripts_path = os.path.join(gencode_tempdir, 'hg38_wgEncodeGencodeBasicV33.txt')
    attributes_path = os.path.join(gencode_tempdir, 'wgEncodeGencodeAttrsV33.txt')
    coding_transcripts_path = os.path.join(gencode_dir, 'hg38_ensGene.txt')
    coding_mrna_path = os.path.join(gencode_dir, 'hg38_ensGeneMrna.fa')

    # 1. Create a file with the coding transcript set
    coding_transcripts = set()
    attributes_df = pd.read_csv(attributes_path, sep='\t', header=None)
    for _, row in attributes_df.iterrows():
        gene_biotype, transcript_id, transcript_biotype, clazz = row[2], row[4], row[6], row[12]
        if gene_biotype == "protein_coding" and transcript_biotype == "protein_coding" and \
                clazz == "coding":
            coding_transcripts.add(transcript_id)

    chromosomes = [str(c) for c in range(1, 23)] + ['X', 'Y']
    chromosomes = ['chr' + c for c in chromosomes]
    transcripts_df = pd.read_csv(transcripts_path, sep='\t', header=None)
    
    with open(coding_transcripts_path, 'w') as coding_transcripts_file:
        for _, row in transcripts_df.iterrows():
            transcript_id, chrom = row[1], row[2]
            if transcript_id in coding_transcripts and chrom in chromosomes:
                line = '\t'.join(map(str, row.tolist()))
                coding_transcripts_file.write(line + '\n')

    # 2. Create a file with mRNA sequences of the coding proteins
    seqfile = os.path.join(seqdir, "hg38.fa")
    bashCmd = "perl {}/retrieve_seq_from_fasta.pl -format genericGene -seqfile {} -outfile {} {}".format(
            annovar_dir, seqfile, coding_mrna_path, coding_transcripts_path)
    subprocess.run(bashCmd.split())

def run(annovar_dir):
    gencode_dir = os.path.join(annovar_dir, "gencode33_ensembl99")
    gencode_tempdir = os.path.join(gencode_dir, "tempdir")
    seqdir = os.path.join(annovar_dir, "hg38", "seq")

    make_directory_if_necessary(gencode_dir)
    make_directory_if_necessary(gencode_tempdir)
    make_directory_if_necessary(seqdir)

    if is_already_prepared(gencode_dir):
        return

    # Preparation involves two steps:
    #   1. Downloading wgEncodeGencodeBasicV33
    #   2. Selecting protein-coding genes
    download_gencode(annovar_dir, seqdir, gencode_tempdir)
    filter_transcript_set(annovar_dir, seqdir, gencode_dir, gencode_tempdir)

if __name__ == '__main__':
    annovar_dir = sys.argv[1]
    run(annovar_dir)
