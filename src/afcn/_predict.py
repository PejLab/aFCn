"""Predict gene expression from phased genotypes.

By: Genomic Data Modeling Lab
"""

import os
import numpy as np

from . import model
from . import bedio
from . import vcfio


FIELD_DELIMITER = "\t"
N_LINES = 500


def run(vcf, par_file, output_dir):

    output_file = os.path.join(output_dir, "predictions.bed")

    reference_expression = 0

    # open all files
    with (bedio.open_param(par_file, "r") as fpars,
          vcfio.read_vcf(vcf) as fvcf,
          bedio.open_predict(output_file, "w") as fout):

        # write meta data to output_file
        fout.meta["vcf"] = vcf
        fout.meta["parameter_file"] = par_file
        fout.sample_names = fvcf.samples

        # Perform gene expression predictions

        for gene_id, variants in fpars.group_by("gene_id"):

            log2_afc = np.zeros(len(variants))
            hap_one = np.zeros(shape=(len(variants), fvcf.n_samples))
            hap_two = np.zeros(shape=(len(variants), fvcf.n_samples))

            # get sample genotypes of each gene associated variant
            for i, v in enumerate(variants):

                log2_afc[i] = v[fpars.idx("log2_afc")]

                tmp = fvcf.get_biallelic_genotypes(v[fpars.idx("chrm")],
                                                   v[fpars.idx("pos")])

                if not tmp["phased"]:
                    raise ValueError(("Data are not phased, "
                                    "as is required."))

                hap_one[i, :] = tmp["genotypes"][0,:]
                hap_two[i, :] = tmp["genotypes"][1,:]

            # given haplotypes and parameters compute expected gene expression
            gene_expr = model.predict(hap_one.T,
                                      hap_two.T,
                                      reference_expression,
                                      log2_afc).astype(str)

            # not all rec values from this iteration of record should
            # have identical genomic coordinates.

            fout.write(rec[fpars.idx("chrm")],
                       v[fpars.idx("gene_start")],
                       v[fpars.idx("gene_end")],
                       gene_id,
                       *gene_expr)
