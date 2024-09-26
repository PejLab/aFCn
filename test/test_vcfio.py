"""Tests of for vcf parsing"""

import os
import unittest
import tempfile

import numpy as np
import pysam

from . import simulate_vcf as simvcf

from afcn import vcfio



def setUpModule():
    global vcf_data 
    global vcf_dir 
    global vcf_name

    vcf_dir = tempfile.TemporaryDirectory()
    vcf_data = simvcf.VCFDataSet()

    tmp_vcf_name = os.path.join(vcf_dir.name, "tmp.vcf")

    with open(tmp_vcf_name, "w") as fid:
        fid.write(str(vcf_data))

    vcf_name = pysam.tabix_index(tmp_vcf_name,
                                 preset="vcf")


def tearDownModule():
    vcf_dir.cleanup()


class TestVCF(unittest.TestCase):
    def setUp(self):
        self.vcf = vcfio.ParseGenotypes(vcf_name, "r")

    def test_properties(self):
        # equal sample size and equal sample names
        self.assertEqual(self.vcf.n_samples,
                         len(vcf_data.samples))

        for test_val,true_val in zip(self.vcf.samples,
                                     vcf_data.samples):
            self.assertEqual(test_val, true_val)

    def test_genotypes_unphased_data(self):
        """Verify genotype records for unphased data."""
        for record in vcf_data.records:

            if record.is_phased():
                continue

            outputs = self.vcf.get_genotypes(record.chrom,
                                         record.pos-1)

            for out in outputs:
                if not record.is_qc_passed():
                    self.assertNotEqual(out["status"], 0)
                    continue

                self.assertFalse(out["phased"])

                self.assertEqual(len(out["alts"]),
                                 len(record.alt))

                for alt_allele in out["alts"]:
                    self.assertIn(alt_allele, record.alt)

                for alt_allele in record.alt:
                    self.assertIn(alt_allele, out["alts"])


                self.assertTupleEqual(record.genotypes.shape,
                                      out["genotypes"].shape)

                for hap_num in range(2):
                    for data_val, true_val in zip(out["genotypes"][hap_num,:],
                                                  record.genotypes[hap_num,:]):
                        if np.isnan(true_val):
                            self.assertTrue(np.isnan(data_val))
                        else:
                            self.assertEqual(data_val, true_val)

    def test_genotype_failed_query(self):
        """Verify that None for filtered variants.
        """
        for record in vcf_data.records:

            outputs = self.vcf.get_genotypes(record.chrom,
                                         record.pos-1)

            for out in outputs:
                if record.alt == '.':
                    self.assertEqual(out["status"], 2)
                elif record.filter != "PASS":
                    self.assertEqual(out["status"], 3)
                else:
                    for allele in out["alts"]:
                        self.assertIn(allele, record.alt)

                    self.assertEqual(out["status"], 0)


    def test_genotype_phased_data(self):
        """Verify genotype records for phased data."""

        for record in vcf_data.records:
            outputs = self.vcf.get_genotypes(record.chrom,
                                         record.pos-1)

            for out in outputs:
                if not record.is_qc_passed():
                    self.assertNotEqual(out["status"], 0)
                    continue

                if not record.is_phased():
                    self.assertFalse(out["phased"])
                    continue

                self.assertTrue(out["phased"])

                for allele in out["alts"]:
                    self.assertIn(allele, record.alt)

                true_genotypes = record.genotypes

                self.assertTupleEqual(out["genotypes"].shape, (2, len(self.vcf.samples)))

                for i in range(len(record.samples)):
                    if np.isnan(true_genotypes[0, i]):
                        self.assertTrue(np.isnan(out["genotypes"][0,i]))
                    else:
                        self.assertEqual(true_genotypes[0, i], out["genotypes"][0, i])

                    if np.isnan(true_genotypes[1, i]):
                        self.assertTrue(np.isnan(out["genotypes"][1,i]))
                    else:
                        self.assertEqual(true_genotypes[1, i], out["genotypes"][1, i])

    def test_no_variants_found(self):
        """Verify genotype records for phased data."""
        # note that in the simulated vcf there is no variant
        # at position 0
    
        for record in vcf_data.records:
            break

        outputs = self.vcf.get_genotypes(record.chrom, 1)
        for out in outputs:
            self.assertEqual(out["status"], 1)


class TestFilters(unittest.TestCase):
    def setUp(self):
        self.vcf = vcfio.ParseGenotypes(vcf_name, "r")
        
    def test_list_requirement(self):
        rec = vcf_data.records[0]

        with self.assertRaises(TypeError):
            next(self.vcf.get_genotypes(rec.chrom,
                                   rec.pos-1,
                                   filter_vals = rec.filter))

    def test_list_inputs(self):
        for rec in vcf_data.records:
            outputs = self.vcf.get_genotypes(rec.chrom,
                                         rec.pos-1,
                                         ['pass','.', 'missing'])

            for out in outputs:
                if ((filt := rec.filter.lower()) == "pass"
                    or filt == "."
                    or filt == "missing"):
                    self.assertEqual(out["status"], 0)
                    continue

                self.assertNotEqual(out["status"], 0)


if __name__ == "__main__":
    unittest.addModuleCleanup(tearDownModule)
    unittest.main()
    unittest.doModuleCleanup()
