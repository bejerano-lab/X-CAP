This folder contains the results of X-CAP scores when evaluated against variants in chromosome 21
of D_validation. X-CAP features can be found at results/stopgains.features, and X-CAP scores
can be found at results/xcap_scores.vcf. 

The following program was used to generate X_CAP scores:
  python run_xcap.py \
      ../example/chr21.vcf  \           # input VCF file
      ../example/ref.vcf    \           # reference VCF file
      ../example/results    \           # output directory
      /absolute/path/to/annovar/dir \   # annovar directory
 
The ref.vcf file has not been included in this zip file because of archive size limitations.

To run X_CAP without having the ref.vcf file, one can run
  python run_xcap.py \
      ../example/chr21.vcf \            # input VCF file
      ../example/empty.vcf \            # reference VCF file
      ../example/results_without_ref \  # output directory
      /absolute/path/to/annovar/dir \   # annovar directory

Note that these results will be worse because X-CAP uses the reference data to make predictions at test
time.
