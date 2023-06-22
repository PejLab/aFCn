import numpy as np
from pysam import VariantFile


# TODO automatic generation of tabix index from bgzip vcf
def read_vcf(vcf):
    """Open file object """
    if not os.exists(f"{vcf}.tbi"):
        raise ValueError("Need tabix index for vcf file.")

    return ParseGenotypes(vcf, "r", index_filename=f"{vcf}.tbi")


# TODO need to figure out indexing rules.
class ParseGenotypes(VariantFile):
    @property
    def samples(self):
        return self.header.samples

    @property
    def n_samples(self):
        return len(self.header.samples)

    def get_biallelic_genotypes(self, contig, pos):
        """Get array(s) that encode genotypes for specified variant.

        For a genotype to be phased, all samples must be phased.

        Args:
            contig: (str) the contig in which the variant is located.
            pos: (int) 

        Returns:
            if not biallelic: None
            dict: {
                    phased: (bool),
                    genotype: ((n_sample,) np.ndarray) if not phased,
                        otherwise ((2, n_sample) np.ndarray)
                   }
        """

        genotypes = np.zeros(shape=(2, self.n_samples))

        phased = True

        for i, variant in enumerate(self.fetch(chrom, pos-1, pos)):

            if len(variant.alts) > 1:
                raise RuntimeWarning(("Not biallelic variant "
                                      f"{chrom}:{start}-{end}, "
                                      "no results computed"))
                return None

            var_encoding = {variant.ref:0,
                            variant.alts[0]:1}

            for n, samp_geno in enumerate(variant.itervalues()):


                genotypes[0, n] = var_encoding[samp_geno.alleles[0]]
                genotypes[1, n] = var_encoding[samp_geno.alleles[1]]

                # if geno.phased is False for any one sample, then by
                # this if statement geno.phased will be False at the
                # termination of for loop
                if phased:
                    phased = samp_geno.phased

        if not phased:
            genotypes = np.sum(genotypes, 0)

        return {"phased":phased,
                "genotypes":genotypes}



