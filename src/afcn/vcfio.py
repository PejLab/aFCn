"""VCF IO tools.

By: GDML
"""
import os
import numpy as np
from pysam import VariantFile


# TODO automatic generation of tabix index from bgzip vcf
def read_vcf(vcf):
    """Open file object """
    if not os.path.exists(f"{vcf}.tbi"):
        raise ValueError("Need tabix index for vcf file.")

    tbx_file = f"{vcf}.tbi"
    if os.path.exists(tbx_file)
        return ParseGenotypes(vcf, "r", index_filename=tbx_file)

    raise RuntimeWarning("No tabix file detected.")
    return ParseGenotypes(vcf, "r")


# TODO need to figure out indexing rules.
class ParseGenotypes(VariantFile):
    @property
    def samples(self):
        return list(self.header.samples)

    @property
    def n_samples(self):
        return len(self.header.samples)

    def get_biallelic_genotypes(self, contig, pos, filter_val="PASS"):
        """Get array(s) that encode genotypes for specified variant.

        For a genotype to be phased, all samples must be phased.

        Args:
            contig: (str) the contig in which the variant is located.
            pos: (int) 

        Returns:
            None: If 
                * variant is not biallelic
                * if FILTER != filter_val
            dict: {
                    phased: (bool),
                    genotype: ((n_sample,) np.ndarray) if not phased,
                        otherwise ((2, n_sample) np.ndarray)
                   }

        Raises:
            RuntimeWarning : If 
        """

        genotypes = np.zeros(shape=(2, self.n_samples))
        phased = True

        for i, variant in enumerate(self.fetch(contig, pos-1, pos)):

            # return none if variant isn't an SNV, and doesn't pass
            # filter, capitilization matters
            if (len(variant.alts) > 1 or
                len(filter_vals := variant.filter.keys()) != 1 or
                filter_vals[0] != filter_val):
                return None

            var_encoding = {variant.ref:0,
                            variant.alts[0]:1}

            for n, samp_geno in enumerate(variant.samples.itervalues()):


                genotypes[0, n] = var_encoding[samp_geno.alleles[0]]
                genotypes[1, n] = var_encoding[samp_geno.alleles[1]]

                # if geno.phased is False for any one sample, then by
                # this if statement geno.phased will be False at the
                # termination of for loop
                if phased:
                    phased = samp_geno.phased

        if not phased:
            genotypes = np.sum(genotypes, 0)
            raise RuntimeWarning(f"Variant {variant.id} not phased.")

        return {"phased":phased,
                "genotypes":genotypes}



