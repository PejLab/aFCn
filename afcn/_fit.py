"""Fit afc model to data

By: Genomic Data Modeling Lab
"""

import os
import logging
import numpy as np

from . import model
from . import bedio
from . import vcfio


def run(vcf, gene_expr, gene_and_variants, output_prefix):

    if output_prefix None:
        output_prefix = "fit"

    output_fname = f"{output_prefix}{bedio.BED_SUFFIX}"

    logging.info("Start: effect size inference")
    logging.info(f"vcf file:{vcf}")
    logging.info(f"expression file:{gene_expr}")
    logging.info(f"variant set per gene:{gene_and_variants}")

    with (vcfio.read_vcf(vcf) as fin_vcf,
          bedio.read_gene_expression(gene_expr) as fin_expr,
          bedio.read_gene_variant_map(gene_and_variants) as fin_gene,
          bedio.open_param(output_fname, "w") as fout_param):
        pass


     

    logging.info(f"output written to:{output_fname}")
    logging.info("Finished")

    raise NotImplementedError
