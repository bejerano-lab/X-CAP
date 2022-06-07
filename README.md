# X-CAP

X-CAP is a GBT model that predicts the pathogenicity of single-nucleotide stopgain variants. This repository contains the trained model as well as code to run X-CAP on a VCF file of variants. 

## Repository structure
- ```bin```: the trained X-CAP model
- ```data```: variants in D_original and D_validation (HGMD variants are only labeled with accession numbers)
- ```example```: demo of X-CAP on chromosome 21 variants in D_validation
- ```figures```: code to generate figures/tables in the paper
- ```predictions```: X-CAP predictions for all stopgains in the human proteome
- ```src```: code to generate X-CAP features and run X-CAP on a VCF file

## Setup
1. Create a Conda environment with necessary requirements.
 ```Python
    conda env create --file environment.yml
 ```
 2. Download [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/user-guide/download/).

## Quick usage guide
To run X-CAP on an arbitrary VCF file, use the following command:
```Python
  python src/run_xcap.py <input_vcf_file> <reference_vcf_file> <output_dir> <annovar_dir>
```
1. ```<input_vcf_file>```: Path to VCF file containing the variants of uncertain significance (VUS) which should be scored. X-CAP will score only those variants annotated as stopgains by ANNOVAR and disregard the rest.
2. ```<reference_vcf_file>```: Path to VCF file containing a reference set of known stopgain variants. X-CAP uses these to produce features for variants in the ```<input_vcf_file>```. During evaluation of X-CAP, we used training variants from D_original as the reference set.
3. ```<output_dir>```: Directory where X-CAP results will be written. X-CAP features will be published to ```<output_dir>/stopgains.features``` and scores to ```<outptut_dir>/xcap_scores.vcf```.
4. ```<annovar_dir>```: Directory containing the ANNOVAR program. X-CAP uses ANNOVAR to filter out non-stopgain variants.

An illustrative demo can be found in the ```examples``` subdirectory.
