"""The afcn program conducts haplotype analysis of gene expression.

By: Genomic Data Modeling Laboratory
"""

LOG_SUFFIX = ".log"

import sys
import os
import argparse
import logging


parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("--version",
                    dest="version",
                    action="store_true")

subparsers = parser.add_subparsers(title = "subcommands",
                                   dest = "subparser_name",
                                   help = "Select the task for afcn to complete.")


# ================================================================
# load input / outoput file specifications

with open(os.path.join(os.path.dirname(__file__),
                       "filespecs", "gene_expression.md"),
          "r") as fid:
    _spec_gene_expression = fid.read()

with open(os.path.join(os.path.dirname(__file__),
                       "filespecs", "param.md"),
          "r") as fid:
    _spec_param = fid.read()

with open(os.path.join(os.path.dirname(__file__),
                       "filespecs", "eqtls.md"),
          "r") as fid:
    _spec_eqtls = fid.read()
          
# ================================================================

# fit subcommand

fit_parser = subparsers.add_parser(
        "fit",
        help="Fit gene expression effect sizes from phased genotypes.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Fit model parameters and write to file in direcotry.",
        epilog =("""

INPUTS
"""
f"\n\n{_spec_gene_expression}\n"
"""
### SPECIFICATION: VCF
    Defined at: https://samtools.github.io/hts-specs/VCFv4.4.pdf
"""
f"\n\n{_spec_eqtls}\n\n"
"""
RETURNS

    <user_prefix>.bed: Parameter bed file
    <user_prefix>.log: error messages
"""
f"\n\n{_spec_param}"))

fit_parser.add_argument("--vcf", 
                        type=str,
                        required=True, 
                        help="Genotype VCF file name")

fit_parser.add_argument("--expr", 
                        type=str,
                        required=True, 
                        help="BED file containing 4 BED fields and "
                            "gene expression per sample.")

fit_parser.add_argument("--eqtls", 
                        type=str,
                        default=None,
                        help=("Name of file that contains eQTL data for which log "
                              "allelic fold change values are inferred from "
                              "data. The first and second columns must be "
                              "#gene_id and variant_id, "
                              "respectively (default None)"))

fit_parser.add_argument("-p",
        type=str,
        default=None,
        help="Prefix of files produced fit submodule.")

fit_parser.add_argument("-o",
                        type=str,
                        default=None,
                        help=("Path and file prefix to print results and logs,"
                              " if specified direcotory(ies) doesn't exist,"
                              " then create."))


# ================================================================

# prediction subcommand



predict_parser = subparsers.add_parser(
        "predict",
        help="Predict gene expression from phased genotype data.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Given an individual's cis-regulatory variant genotype
and log2 fold change effect sizes predict gene expression,
in log2 scale, under the model by Mohammadi et al.
Genome Research 2017.""",
        epilog =("""

RETURNS

    <user_prefix>.bed: expression bed file containing gene expression
        predictions per gene per sample
    <user_prefix>.log: error messages

default prefix *user_prefix* = 'predict'

### SPECIFICATION: vcf
    Defined at: https://samtools.github.io/hts-specs/VCFv4.4.pdf
"""
f"\n\n{_spec_param}\n\n{_spec_gene_expression}"))

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
predict_parser.add_argument("--filters",
            default=["PASS"],
            type=str,
            nargs="+",
            help="""Filter(s) that genotypes must meet to be
                 included in the analysis.  Note that if you 
                 would like to include variants for which the VCF
                 filter value is '.', include 'missing' as a
                 filter value.  Input is not
                 case sensitive, default pass""")

predict_parser.add_argument("-p",
        type=str,
        default=None,
        help="Prefix of files produced fit submodule.")

predict_parser.add_argument(
        "-o",
        type=str,
        default=None,
        help="Directory to print results, if does not exists, create.")

predict_parser.add_argument(
        "--scale",
        type=str,
        nargs="?",
        default="linear",
        choices= ["linear", "log2", "log", "log10"],
        help="Report predictions in linear, log2, log10, or log (natural log) scale")

predict_parser.add_argument(
        "--decimals",
        type=int,
        default=None,
        help=("Number of decimal places to round predictions."
              " Rounding is performed in the scale of the returned data."))



# ================================================================

# parse input

args = parser.parse_args(sys.argv[1:])

if args.version:
    import afcn
    print(f"afcn {afcn.__version__}")
    quit()


# default directories and prefix for output files
if args.o is None:
    args.o = f"afcn_{args.subparser_name}"

if args.p is None:
    args.p = "output"

# remove and path related information
args.p = os.path.basename(args.p)


# make all directories that do not exist 
out_path = None

for directory in args.o.strip(os.path.sep).split(os.path.sep):
    if out_path is None:

        out_path = directory

    else:

        out_path = os.path.join(out_path, directory)

    if not os.path.exists(out_path):
        os.mkdir(out_path)


args.o = os.path.join(out_path, args.p)

if os.path.sep in args.o and not os.path.exists(os.path.dirname(args.o)):
    os.mkdir(os.path.dirname(args.o))

logging.basicConfig(filename=f"{args.o}{LOG_SUFFIX}",
                level=logging.INFO,
                filemode="w",
                format="%(levelname)s\t%(asctime)s\t%(message)s",
                datefmt="%Y-%m-%d %H:%M:%S")
logging.info(" ".join(sys.argv))


if args.subparser_name == "fit":
    from . import _fit

    _fit.run(args.vcf, args.expr, args.eqtls, args.o)


if args.subparser_name == "predict":
    from . import _predict
    _predict.run(args.vcf, args.params, args.o,
                 args.filters, args.scale, args.decimals)

