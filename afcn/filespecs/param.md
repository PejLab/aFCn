### SPECIFICATION: PARAMETER FILE

BED 4 file with meta data and header as follows:

* line [0, n): meta data prepended with ##
* line n: header prepended with #
* line [n+1, n+1+N_records): records

Contents:

* meta data
    - ##afcn_version = (version str)
    - ##vcf_file = (str)
    - ##expression_file = (str)
    - ##eqtl_file = (str)
    - etc.
* BED fields labeled: 
    - #chrom, 
    - qtl_start: (int) minimum(variant position, gene position)
    - qtl_end: (int) maximum(variant position, gene position)
    - qtl_id: (string) gene_id/variant_id
* custom fields:
    - gene_start: (integer) start of genomic feature
    - gene_end: (integer) end of genomic feature
    - gene_id: (string) gene name
    - variant_pos: (integer) variant genomic coordinates from VCF
    - variant_id: (string)
    - ref: (char) reference allele
    - alt: (char) alternative allele
    - log2_afc: (float) estimated model parameter
    - remaining columns are statistics relevant to
        parameter inference
