"""Predict gene expression from phased genotypes.

By: Genomic Data Modeling Lab
"""

import os
import logging
import numpy as np

from . import model
from . import dataio


FIELD_DELIMITER = "\t"
N_LINES = 500

logging.basicConfig(filename=os.path.join(args.output_dir, 
                                          "predictions.log"),
                    format="%(asctime)s %(levelname)s : %(message)s",
                    level = logging.INFO)



def run(vcf, par_file, output_dir):

    logging.INFO("Start: gene expression predictions")
    logging.INFO(f"vcf file:{vcf}")
    logging.INFO(f"parameter file:{par_file}")

    output_file = os.path.join(output_dir, "predictions.bed")

    # open all files
    with (dataio.read_bed(par_file) as fpars,
          dataio.read_vcf(vcf) as fvcf,
          dataio.open_bed(output_file, "w") as fout):

        #vcf_pars = dataio.GenotypeVcfParser(fvcf)

        for gene_id, records in fpars.group_by("gene_id"):

            log2_afc = np.zeros(len(records))

            #for i, record in enumerate(records):
            #    log2_afc[i] = record[fpars.idx("log2_afc")]
            #    
            #    hap_one, hap_two = fvcf.fetch(record[fpars.idx("chrm")],
            #                            record[fpars.idx("pos")])



            # compute gene expression






    # organize data to be compatible with model.predict
    # make prediction, add to buffer

    logging.INFO(f"output written to:{output_dir/predictions.bed}")
    logging.INFO("Finished")


