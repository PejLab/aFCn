"""The aFCn program conducts haplotype analysis of gene expression.

By: Genomic Data Modeling Laboratory
"""

import sys
import argparse


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("--version",
                    dest="version",
                    action="store_true")

subparsers = parser.add_subparsers(title = "subcommands",
                                   dest = "subparser_name",
                                   help = "Select the task for afcn to complete.")


# ================================================================

# fit subcommand

# fit_parser = subparsers.add_parser(
#         "fit",
#         help="Fit gene expression effect sizes from phased genotypes.",
#         formatter_class=argparse.RawDescriptionHelpFormatter,
#         description="Fit model parameters and write to file in direcotry.",
#         epilog = """OUTPUT FILE SPECIFICATION
# """)
# 
# fit_parser.add_argument("--vcf", 
#                         type=str,
#                         required=True, 
#                         help="Genotype VCF file name")
# 
# fit_parser.add_argument("--expr", 
#                         type=str,
#                         required=True, 
#                         help="BED file containing 4 BED fields and "
#                             "gene expression per sample.")
# 
# 
# # optional inputs
# 
# fit_parser.add_argument("--eqtl", 
#                         type=str,
#                         default=None,
#                         help=("Name of file that contains eQTLs for which log "
#                               "allelic fold change values are inferred from "
#                               "data. The first and second columns must be "
#                               "#gene_id and variant_id, "
#                               "respectively (default None)"))
# 
# 
# # TODO write specification of covariate file
# 
# fit_parser.add_argument("--covariates",
#                         type=str,
#                         required=False,
#                         help=("File containing covariate data."))
# 
# fit_parser.add_argument("output_dir",
#                         type=str,
#                         nargs="?",
#                         default="afcs.bed", 
#                         help="Output file name")
# 

# ================================================================

# prediction subcommand

predict_parser = subparsers.add_parser(
        "predict",
        help="Predict gene expression from phased genotype data.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Given an individual's cis-regulatory variant genotype
and log2 fold change effect sizes predict gene expression 
under the model by Mohammadi et al. Genome Research 2017.""",
        epilog = """
SPECIFICATION PARAMETER FILE INPUT

BED 4 file with meta data and header as follows:

* line [0, n): meta data prepended with ##
* line n: header prepended with #
* line [n+1, n+1+N_eqtls): records

Contents:

## meta data
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


OUTPUT FILES

* <output_dir>/predictions.bed
* <output_dir>/predictions.log


SPECIFICATION predictions.bed

BED 4 file with meta data and header as follows:

* line [0, n): meta data prepended with ##
* line n: header prepended with #
* line [n+1, n+1+N_samples): records

Contents:

* meta data
    - ##afcn_version=(version str)
    - ##vcf_file=(str)
    - ##parameter_file=(str)
* BED fields labeled: 
    - #chrom, start, end, name
* custom fields:
    - sample_id: (float) predicted gene expression per sample
""")

predict_parser.add_argument(
        "--vcf",
        required=True,
        help="""VCF file version 4.4 as defined in the samtools 
        documentation.""")
predict_parser.add_argument(
        "--params",
       required=True,
       help="""Tab delimited text file providing
       log2 aFC point estimates per (gene, variant)
       pair.  See below for more details.
       """)
predict_parser.add_argument(
        "output_dir",
        type=str,
        help="Directory to write results")

# ================================================================

# twas subcommand

# twas_parser = subparsers.add_parser(
#         "twas",
#         help="""Perform transcriptome wide association study
#         from gene expression predictions.""",
#         formatter_class=argparse.RawDescriptionHelpFormatter,
#         description="""
# Given a set of gene expression predictions per individual, perform 
# association tests between the predicted expression and a observed
# phenotype.
# """,
#         epilog = """
# OUTPUT FILE SPECIFICATION
# """)


args = parser.parse_args(sys.argv[1:])


if args.version:
    import afcn
    print(f"afcn {afcn.__version__}")

if args.subparser_name == "fit":
    from . import _fit

    _fit.run()


if args.subparser_name == "predict":
    from . import _predict

    _predict.run(args.vcf, args.params, args.output_dir)


if args.subparser_name == "twas":
    from . import _twas

    _twas.run()
