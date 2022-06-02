import argparse
import lightgbm as lgb
import os
import sys

import prepare_annovar
import filter_stopgains

sys.path.append(os.path.join(sys.path[0], "utils"))
sys.path.append(os.path.join(sys.path[0], "featurization"))
import featurize
import evaluation_utils
import annovar_utils

XCAP_THRESHOLD = 0.06011546409371992

def parse_arguments():
    parser = argparse.ArgumentParser(description="runs X-CAP on an arbitrary VCF file")
    parser.add_argument("input_vcf_file", type=str, help="input VCF file")
    parser.add_argument("reference_vcf_file", type=str, help="training set VCF file")
    parser.add_argument("output_dir", type=str)
    parser.add_argument("annovar_dir", type=str)
    parser.add_argument("--classifier_path", type=str, help="path to .mdl file",
            default=os.path.join(sys.path[0], "..", "bin", "xcap_052021.mdl"))
    args = parser.parse_args()
    return args

def write_vcf_header(output_f):
    output_f.write("##fileformat=VCFv4.1\n")
    output_f.write("##reference=hg38\n")
    columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    output_f.write('#' + '\t'.join(columns) + '\n')

def get_predicted_label(score):
    return int(score > XCAP_THRESHOLD)

def publish_results(output_path, stopgains_path, features_path, predictions):
    feature_names, features, _ = evaluation_utils.get_xcap_data(features_path)

    with open(output_path, 'w') as output_f:
        write_vcf_header(output_f)

        for i, line in enumerate(open(stopgains_path)):
            v = annovar_utils.get_variant_from_line(line)
            
            info = line.strip().split('\t')[-1]
            info += ";SCORE={:.4f}".format(predictions[i])
            
            predicted_label = get_predicted_label(predictions[i])
            info += ";PREDICTED_LABEL={}".format(predicted_label)
            
            row = ['chr' + v.chrom, v.pos, '.', v.ref, v.alt, '.', '.', info]
            output_f.write('\t'.join(map(str, row)) + '\n')

def main():
    args = parse_arguments()
    
    # Prepare ANNOVAR to support protein-coding GencodeV33 transcript set
    prepare_annovar.run(args.annovar_dir)

    # Use ANNOVAR to extract stopgains from the VCF files of both the input and the reference.
    # The training set, by construction, only contains stopgains so this procedure is only useful
    # for converting the VCF file into a format that can be parsed by the featurizer.
    input_annovar_annotations_dir = os.path.join(args.output_dir, "input_annovar_files")
    ref_annovar_annotations_dir = os.path.join(args.output_dir, "ref_annovar_files")

    input_stopgains_path = filter_stopgains.run(
            args.input_vcf_file, input_annovar_annotations_dir, args.annovar_dir)
    ref_stopgains_path = filter_stopgains.run(
            args.reference_vcf_file, ref_annovar_annotations_dir, args.annovar_dir)

    # Create features for each stopgain in stopgains_path
    features_path = os.path.join(args.output_dir, "stopgains.features")
    featurize.run_on_file(
            input_stopgains_path, ref_stopgains_path, features_path, args.annovar_dir, False)

    # Get X-CAP predictions and corresponding percentiles
    clf = lgb.Booster(model_file=args.classifier_path)
    _, predictions = evaluation_utils.get_xcap_predictions(clf, features_path)

    # Publish results to another VCF file
    output_path = os.path.join(args.output_dir, "xcap_scores.vcf")
    publish_results(output_path, input_stopgains_path, features_path, predictions)

if __name__ == '__main__':
    main()
