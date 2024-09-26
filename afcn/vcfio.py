"""VCF IO tools.

By: GDML
"""
import os
import numpy as np
from pysam import VariantFile



# TODO automatic generation of tabix index from bgzip vcf
def read_vcf(vcf):
    """Open file object """
    tbx_file = f"{vcf}.tbi"
    if not os.path.exists(tbx_file):
        raise ValueError("Need tabix index for vcf file.")

    if os.path.exists(tbx_file):
        return ParseGenotypes(vcf, "r", index_filename=tbx_file)

    raise RuntimeWarning("No tabix file detected.")


# TODO confirm indexing rules.
class ParseGenotypes(VariantFile):
    @property
    def samples(self):
        return list(self.header.samples)

    @property
    def n_samples(self):
        return len(self.header.samples)

    @property
    def contigs(self):
        return list(self.header.contigs)

    def len_contig(self, contig):
        return self.header.contigs.get(contig).length

    def get_contig_variants(self, contig):
        positions = dict()

        for record in self.fetch(contig):
            positions[record.id] = record.pos

        return positions

    def _extract_genotype_from_pysam_record(self, record, filter_criteria):
        # return none if variant isn't an SNV, and doesn't pass
        # filter, capitilization matters:
        #   * deletion, e.g. Ref field = "."
        #   * more than one filter value
        #   * specified filter satsified

        genotypes = np.full((2, self.n_samples), np.nan)

        if record.alts is None:
            return dict(status=2,
                        msg=("NoAltAllele"
                             f" VCF locus {record.contig}:{record.pos}"))

        # the lower operation ensures that pattern matching
        # is not case sensitive
        satisfy_filter_criterion = False
        record_filter_vals = []
        for w in record.filter.keys():
            record_filter_vals.append(w.lower())

        if len(record_filter_vals) == 0 :
            record_filter_vals.append("missing")

        for fval in filter_criteria:
            if fval.lower() in record_filter_vals:
                satisfy_filter_criterion = True
                break

        if not satisfy_filter_criterion:
            return dict(status=3,
                        msg=("FilterMismatch VCF locus"
                             f" {record.contig}:{record.pos}"
                             f" filter values ({','.join(record_filter_vals)}),"
                             f" valid filter values ({','.join(filter_criteria)})"))


        # according to vcf4.4 specification
        genotype_map = {record.ref:0,
                        None:np.nan}


        for alt_allele in record.alts:
            genotype_map[alt_allele] = None

        phased = True

        for n, samp_geno in enumerate(record.samples.itervalues()):

            # sets the alternative allele encodings
            for allele, idx in zip(samp_geno.alleles,
                                   samp_geno.allele_indices):


                # default entry is np.nan
                if allele is None and idx is None:
                    continue

                elif genotype_map[allele] is None:
                    genotype_map[allele] = idx

                elif genotype_map[allele] != idx:
                    raise ValueError("Inconsistent indexing of"
                                     " alleles across samples")


            genotypes[0, n] = genotype_map[samp_geno.alleles[0]]
            genotypes[1, n] = genotype_map[samp_geno.alleles[1]]

            # if geno.phased is False for any one sample, then by
            # this if statement geno.phased will be False at the
            # termination of for loop
            if phased:
                phased = samp_geno.phased

        return {"phased":phased,
                "genotypes":genotypes,
                "genotype_map":genotype_map,
                "ref":record.ref,
                "alts":record.alts,
                "status": 0,
                "msg": "success"}


    def get_genotypes(self, contig, pos, filter_vals=["PASS"]):
        """Get array(s) of genotypes.

        Extract the genotype{0,1,2,...,np.nan} for all samples
        of at a specified locus. For a genotypes of a variant
        to be considered phased, all samples must be phased.

        Args:
            contig: (str)
                the contig in which the variant is located.
            pos: (int)
                assume 0 based indexing
            filter_vals:(list)
                list of strings specifing acceptable
                variant filter value(s)

        Returns:
            list, one item per variant found at given locus,
                of dictionaries as specified below:

                dict: 
                    {
                        phased: (bool),
                        genotypes: ((2, n_sample) np.ndarray)
                        var_encodings: (dict)
                            key value pairs of allele (A,T,C,G)
                            and numeric value.  Missing genotypes
                            are None
                        status: (int)
                        msg: (str)
                    }

        where `status` and `msg` reports whether a variant that
        satisfies our requirements is found, if not it reports
        the incompatibility.

        | status code    |  msg class                     |
        | -------------- | ------------------------------ |
        | 0              | Success                        |
        | 1              | VariantNotFound                |
        | 2              | NoAltAllele                    |
        | 3              | FilterMismatch                 |
        | 4              | ValueError in fetch            |
        """

        if not isinstance(filter_vals, list):
            raise TypeError("Input filter value must be a Python list")

        
        # Set i to None to detect when no variant is found.  When no
        # variant is found the for loop below does not assign a value
        # to i.
        msg = dict(status=1,
                    msg=("VariantNotFound at locus"
                    f" {contig}:{pos} not in VCF."))
        try:

            records = self.fetch(contig, pos, pos+1)

        except ValueError as err:

            records = []
            msg = dict(status=4, msg=repr(err))

        rec = None

        for rec in records:

            yield self._extract_genotype_from_pysam_record(rec, filter_vals)

        if rec is None:

            yield msg
