import os
import subprocess
import sys

def make_directory_if_necessary(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

# Pre-condition: prepare_annovar.run() has already been called
def run(input_vcf_file, output_dir, annovar_dir):
    make_directory_if_necessary(output_dir)

    # 1. Convert VCF to avinput file
    avinput = os.path.join(output_dir, "avinput")
    bashCmd = "perl {}/convert2annovar.pl -format vcf4 {} --includeinfo --outfile {}".format(
            annovar_dir, input_vcf_file, avinput)
    subprocess.run(bashCmd.split())

    # 2. Annotate each variant
    annotation_prefix = "annovar"
    bashCmd = "perl {}/annotate_variation.pl -out {}/{} -build hg38 -dbtype ensGene -precedence splicing,exonic {} {}/gencode33_ensembl99".format(
            annovar_dir, output_dir, annotation_prefix, avinput, annovar_dir)
    subprocess.run(bashCmd.split())

    # 3. Select stopgains
    exonic_variants_path = "{}/{}.exonic_variant_function".format(output_dir, annotation_prefix)
    stopgains_path = os.path.join(output_dir, "stopgains.tsv")
    bashCmd = 'grep stopgain {}'.format(exonic_variants_path)

    with open(stopgains_path, 'w') as f:
        child = subprocess.Popen(bashCmd.split(), stdout=f)
        child.wait()

    return stopgains_path

if __name__ == '__main__':
    input_vcf_file = sys.argv[1]
    output_dir = sys.argv[2]
    annovar_dir = sys.argv[3]

    run(input_vcf_file, output_dir, annovar_dir)
