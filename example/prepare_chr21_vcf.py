import sys
sys.path.append("../src/utils")
import annovar_utils

D_validation_input_path = "/cluster/u/rrastogi/ECNN/results/independent_dataset/preprocess/remove_train/test.tsv"
output_path = "chr21.vcf"

def write_vcf_headers(output_file):
    output_file.write("##fileformat=VCFv4.1\n")
    output_file.write("##reference=hg38\n")
    columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    output_file.write('#' + '\t'.join(columns) + '\n')

def main():
    input_file = open(D_validation_input_path)
    output_file = open(output_path, 'w')

    write_vcf_headers(output_file)

    for line in input_file:
        row = line.strip().split('\t')
        
        chrom, pos, ref, alt = row[3], row[4], row[6], row[7]
        
        if not chrom.startswith('chr'):
            chrom = 'chr' + chrom
        if chrom != 'chr21':
            continue
        
        label = annovar_utils.get_value_from_info(row[-1], "LABEL")
        info = "LABEL={}".format(label)
        if label == "0":
            # if benign, also add NHOMALT and AC values
            nhomalt = annovar_utils.get_value_from_info(row[-1], "NHOMALT")
            info += ";NHOMALT={}".format(nhomalt)

            ac = annovar_utils.get_value_from_info(row[-1], "AC")
            info += ";AC={}".format(ac)
        
        output = [chrom, pos, '.', ref, alt, '.', '.', info]
        output_file.write('\t'.join(output) + '\n')
    
    output_file.close()

if __name__ == '__main__':
    main()
