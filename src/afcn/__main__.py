"""The aFCn program conducts haplotype analysis of gene expression.

By: Genomic Data Modeling Laboratory
"""


_fit_descpt="""Fit log fold changes from phased genotypes and expression data.

Output File Specification:
    * tab delimited BED file as defined in

        https://github.com/samtools/hts-specs/

        repository at commit 9ddbc52
    * Contains header
    * Data contains 4 BED fields, and custom fields.  The genomic position is
        that of variant:

        #chrom  start   end variant_id  gene_id ln_afc  ci_low  ci_hi   pval
"""

import sys
import argparse

# from . import model


parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers(title = "subcommands",
                                   dest = "subparser_name",
                                   help = "Select the task for afcn to complete.")

# fit subcommand options
# required inputs

fit_parser = subparsers.add_parser("fit",
                                   help=_fit_descpt)

fit_parser.add_argument("--vcf", 
                        type=str,
                        required=True, 
                        help="Genotype VCF file name")

fit_parser.add_argument("--expr", 
                        type=str,
                        required=True, 
                        help="BED file containing 4 BED fields and "
                            "gene expression per sample.")


# optional inputs

fit_parser.add_argument("--eqtl", 
                        type=str,
                        default=None,
                        help=("Name of file that contains eQTLs for which log "
                              "allelic fold change values are inferred from "
                              "data. The first and second columns must be "
                              "#gene_id and variant_id, "
                              "respectively (default None)"))


# TODO write specification of covariate file

fit_parser.add_argument("--covariates",
                        type=str,
                        required=False,
                        help=("File containing covariate data."))

fit_parser.add_argument("output_file",
                        type=str,
                        nargs="?",
                        default="afcs.bed", 
                        help="Output file name")

# predict subcommand options

predict_parser = subparsers.add_parser("predict",
                                       help="test")


# twas subcommand options
twas_parser = subparsers.add_parser("twas",
                                    help="test")



args = parser.parse_args(sys.argv[1:])

if args.subparser_name == "fit":
    from . import _fit

    _fit.run()


if args.subparser_name == "predict":
    from . import _predict

    _predict.run()


if args.subparser_name == "twas":
    from . import _twas

    _twas.run()
