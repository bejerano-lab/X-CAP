import argparse
import os
import sys
from tqdm import tqdm

import conservation
import gene_essentiality
import locations
import metadata
import omim
import read_through
import reinitiation
import rvis
import splicing
import transcript
import transcript_class
import utils

def init(reference_data, annovar_dir):
    gencode_dir = os.path.join(annovar_dir, "gencode33_ensembl99")
    gencode_transcripts_path = os.path.join(gencode_dir, "hg38_ensGene.txt")
    gencode_mrna_path = os.path.join(gencode_dir, "hg38_ensGeneMrna.fa")
    gencode_attributes_path = os.path.join(gencode_dir, "tempdir", "wgEncodeGencodeAttrsV33.txt")

    utils.init_transcript_map(gencode_transcripts_path, gencode_mrna_path)
    gene_essentiality.init(reference_data)
    rvis.init()
    omim.init()
    transcript_class.init(reference_data)
    metadata.init_variant_data(reference_data)
    splicing.init(gencode_attributes_path)
    conservation.init_bws()

def run_on_data(input_data, reference_data, annovar_dir, is_training_set):
    init(reference_data, annovar_dir)
    header, feature_table, labels = None, [], []
   
    print("Featurizing data")
    for row in tqdm(input_data):
        features, names = [], []

        #### GENE ESSENTIALITY
        gene_oe = gene_essentiality.get_average_gene_oe_lof(row)
        transcript_oe = gene_essentiality.get_average_transcript_oe_lof(row)
        rvis_percentile = rvis.get_average_rvis(row)
        recessive_gene = omim.on_recessive_gene(row)
        dominant_gene = omim.on_dominant_gene(row)
        monoclass_pathogenic_transcript = transcript_class.on_monoclass_pathogenic_transcript(row, is_training_set)
        monoclass_pathogenic_exon = transcript_class.on_monoclass_pathogenic_exon(row, is_training_set)

        features += [gene_oe, transcript_oe, rvis_percentile, recessive_gene, dominant_gene,
                     monoclass_pathogenic_transcript, monoclass_pathogenic_exon]
        names += ["gene_oe_lof", "transcript_oe_lof", "rvis", "recessive_gene", "dominant_gene",
                  "monoclass_pathogenic_transcript", "monoclass_pathogenic_exon"]

        #### ZYGOSITY
        zygosity = metadata.get_zygosity(row)
        features += [zygosity]
        names += ["zygosity"]

        #### VARIANT LOCATION
        distance_from_cds_start = locations.get_distance_from_cds_start(row)
        distance_from_cds_end = locations.get_distance_from_cds_end(row)
        relative_location_in_cds = locations.get_relative_location_in_cds(row)
        features += [distance_from_cds_start, distance_from_cds_end, relative_location_in_cds]
        names += ["distance_from_cds_start", "distance_from_cds_end", "relative_location_in_cds"]

        distance_from_exon_start = locations.get_distance_from_exon_start(row)
        distance_from_exon_end = locations.get_distance_from_exon_end(row)
        relative_location_in_exon = locations.get_relative_location_in_exon(row)
        exon_length = locations.get_overlapped_exon_length(row)
        features += [distance_from_exon_start, distance_from_exon_end, relative_location_in_exon, exon_length]
        names += ["distance_from_exon_start", "distance_from_exon_end", "relative_location_in_exon", "exon_length"]

        overlapped_exon_num = locations.get_overlapped_exon_num(row)
        num_exons = transcript.get_num_exons(row)

        on_autosome = locations.on_autosomal_chromosome(row)
        on_X = locations.on_X_chromosome(row)
        on_Y = locations.on_Y_chromosome(row)

        features += [overlapped_exon_num, num_exons, on_autosome, on_X, on_Y]
        names += ["overlapped_exon_num", "num_exons", "on_autosome", "on_X", "on_Y"]

        #### NMD
        dist_from_last_exon_exon_junction = locations.get_distance_from_last_exon_exon_junction(row)
        percentage_of_transcripts_with_nmd = locations.get_percentage_of_transcripts_with_nmd(row)

        features += [dist_from_last_exon_exon_junction, percentage_of_transcripts_with_nmd]
        names += ["dist_from_last_exon_exon_junction", "%_transcripts_with_nmd"]

        #### SPLICING
        can_be_spliced_out = splicing.can_be_spliced_out(row)
        features += [can_be_spliced_out]
        names += ["can_be_spliced_out"]

        #### STOP CODON READ THROUGH
        kmers, stop_codon_encoding = read_through.get_stop_codon(row, num_context=0)
        features += stop_codon_encoding
        names += kmers

        #### TRANSLATION REINITIATION
        dist_to_next_start_codon = reinitiation.get_dist_to_next_start_codon(row)
        features += [dist_to_next_start_codon]
        names += ["dist_to_next_start_codon"]

        #### CONSERVATION
        overlapped_exon_phylop = conservation.get_mean_phylop_on_overlapped_exon(row)
        overlapped_exon_phastcons = conservation.get_mean_phastcons_on_overlapped_exon(row)
        upstream_phylop = conservation.get_mean_phylop(row, upstream=True)
        upstream_phastcons = conservation.get_mean_phastcons(row, upstream=True)
        downstream_phylop = conservation.get_mean_phylop(row, upstream=False)
        downstream_phastcons = conservation.get_mean_phastcons(row, upstream=False)

        features += [overlapped_exon_phylop, overlapped_exon_phastcons, upstream_phylop, upstream_phastcons,
                     downstream_phylop, downstream_phastcons]
        names += ["overlapped_exon_phylop", "overlapped_exon_phastcons", "upstream_phylop",
                  "upstream_phastcons", "downstream_phylop", "downstream_phastcons"]

        header = names
        feature_table.append(features)
        labels.append(metadata.get_label(row))

    return header, feature_table, labels

def run_on_file(input_filepath, reference_filepath, output_filepath, annovar_dir, is_training_set):
    input_data = [line.strip().split('\t') for line in open(input_filepath)]
    reference_data = [line.strip().split('\t') for line in open(reference_filepath)]
    
    header, features, labels = run_on_data(input_data, reference_data, annovar_dir, is_training_set)
    with open(output_filepath, 'w') as outfile:
        header.append("label")
        outfile.write('\t'.join(header) + '\n')
        for f, l in zip(features, labels):
            output = f + [l]
            output = '\t'.join(map(str, output))
            outfile.write(output + '\n')

def main():
    parser = argparse.ArgumentParser(description="X-CAP feature generation tool")
    parser.add_argument("input_file", type=str, help="file containing the data to featurize")
    parser.add_argument("reference_file", type=str, help="file containing the training dataset")
    parser.add_argument("output_file", type=str, help="file where features will be output")
    parser.add_argument("annovar_dir", type=str, help="ANNOVAR dir prepared with GencodeV33 transcript set")
    parser.add_argument("--is_training_set", type=bool, nargs='?', default=False, const=True,
                       help="input_filepath contains the training_data")
    args = parser.parse_args()

    print("Running X-CAP. Input: {}, reference: {}, output: {}, is training set?: {}".format(
        args.input_file, args.reference_file, args.output_file, args.is_training_set))
    run_on_file(args.input_file, args.reference_file, args.output_file, args.annovar_dir, args.is_training_set)

if __name__ == '__main__':
    main()
