### SPECIFICATION: GENE EXPRESSION BED

BED 4 file with meta data and header as follows:

* line [0, n): meta data prepended with ##
* line n: header prepended with #
* line [n+1, n+1+N_samples): records

Contents:

* Required meta data
    - ##afcn_version=(version str)
    - ##vcf_file=(str)
    - ##parameter_file=(str)
    - ##prediction_scale=(str) whether the data is in linear scale
        or transformed to log (natural log), log2, log10 scale
    - ##prediction_values=(str) describes what type and format of
        data is written in file
    - ##rounding=(str) whether the data were rounded, and if
        so to what decimal place
* BED fields labeled: 
    - chrom
    - start : gene start position
    - end : gene end position
    - name : ensembl gene id
* custom fields:
    - sample_id: (float) per haplotype predicted gene expression,
        linear scale.  e.g. 322|32
