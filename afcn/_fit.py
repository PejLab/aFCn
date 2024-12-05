"""Fit afc model to data

By: Genomic Data Modeling Lab
"""

import os
import logging
import numpy as np

from . import model
from . import bedio
from . import vcfio


def run(vcf, gene_expr_file, eqtl_file, output_prefix):

    if output_prefix None:
        output_prefix = "fit"

    output_fname = f"{output_prefix}{bedio.BED_SUFFIX}"

    logging.info("Start: effect size inference")
    logging.info(f"vcf file:{vcf}")
    logging.info(f"expression file:{gene_expr_file}")
    logging.info(f"eqtl file:{eqtl_file}")

    with (vcfio.read_vcf(vcf) as fin_vcf,
          bedio.read_gene_expression(gene_expr_file) as fin_expr,
          bedio.read_eqtl_map(eqtl_file) as fin_eqtl,
          bedio.open_param(output_fname, "w") as fout_param):

        # meta data
        fout_param.meta["vcf"] = vcf
        fout_param.meta["eqtl_file"] = eqtl_file
        fout_param.meta["gene_expr_file"] = gene_expr_file


        fout_param.write_meta_and_header()

        # loop over all eqtls for a single gene, use gene_id
        for gene_id, variants in fin_eqtl.group_by("gene_id"):

            haplotypes = [np.full((fvcf.n_samples, len(variants)), np.nan),
                          np.full((fvcf.n_samples, len(variants)), np.nan)]

            # get sample genotypes of each gene associated variant
            for i, v in enumerate(variants):
                sample_genotype_records = fvcf.get_genotypes(
                                                v[fpars.idx("chrom")],
                                                v[fpars.idx("variant_pos")],
                                                filter_vals=filters)
            fin_expr

     
            fit_out = model.fit(haplotypes[0], haplotypes[1],
                                gene_expr,
                                reg = reg,
                                reg_const = reg_const)



    logging.info(f"output written to:{output_fname}")
    logging.info("Finished")

    raise NotImplementedError
