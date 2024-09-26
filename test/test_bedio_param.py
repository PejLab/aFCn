import os
import shutil
import io
import tempfile
import time

import unittest

from . import simulate_bed as sim_bed

from afcn import bedio


def setUpModule():
    global bedio_dir
    global bedio_data
    global bedio_gzdata
    global no_header_meta
    global wrong_header_order

    bedio_dir = tempfile.TemporaryDirectory()

    bedio_data = sim_bed.ParamBedDataSet(os.path.join(bedio_dir.name,
                                                      "tmp_param.bed"))
    bedio_gzdata = sim_bed.ParamBedDataSet(os.path.join(bedio_dir.name,
                                                      "tmp_param.bed.gz"),
                                           compress=True)

    no_header_meta = sim_bed.ParamBedDataSet(os.path.join(bedio_dir.name,
                                                "no_header_meta_entries.bed"),
                                        meta=None, header=None)

    wrong_header_order = sim_bed.ParamBedDataSet(os.path.join(bedio_dir.name,
                                                    "wrong_entries.bed"),
                                                 header=[sim_bed.DEFAULT_HEADER[1],
                                                         sim_bed.DEFAULT_HEADER[0],
                                                         *sim_bed.DEFAULT_HEADER[2:]])


def tearDownModule():
    bedio_dir.cleanup()


class Testopen_param(unittest.TestCase):

    def test_correct_inputs(self):
        for ds in (bedio_data, bedio_gzdata):

            with bedio.open_param(ds.filename, "r") as tmp:
                self.assertIsInstance(tmp, bedio.ParseParamBed)
    
    def test_not_bed_exception(self):
        
        with self.assertRaises(ValueError):
            bedio.open_param(os.path.join(bedio_dir.name,
                                          "test.bed.bz2"),
                             "r")

        with self.assertRaises(ValueError):
            bedio.open_param(os.path.join(bedio_dir.name,
                                          "test.vcf.gz"), "r")

        with self.assertRaises(ValueError):
            bedio.open_param(os.path.join(bedio_dir.name, "bed.txt"), "r")

        with self.assertRaises(ValueError):
            bedio.open_param(os.path.join(bedio_dir.name, "test_bed.gz"), "r")


class TestParameterFileParser(unittest.TestCase):

    def test_file_build_fileio_data(self):
    
        for ds in (bedio_data, bedio_gzdata):
    
            with bedio.open_param(ds.filename, "r") as fpar:
    
                for key, meta_val in ds.meta.items():
                    self.assertEqual(meta_val, fpar.meta[key])
    
                for i, header_val in enumerate(ds.header):
                    self.assertEqual(header_val, fpar.header[i])
    
    def test_group_by_fileio_data(self):
    
        for ds in (bedio_data, bedio_gzdata):
    
            with bedio.open_param(ds.filename, "r") as fpar:
    
                feature_idx = 0
                record_idx = 0
    
                for feature, feature_data in fpar.group_by("gene_id"):
                    self.assertEqual(feature, ds.features[feature_idx])
                    feature_idx += 1
    
                    for var_feature_data in feature_data:
    
                        for true_val, parsed_val in zip(ds.data[record_idx],
                                                        var_feature_data):
    
                            self.assertEqual(true_val, parsed_val)
    
                        record_idx += 1
    
    def test_invalid_headers(self):
    
        for ds in (wrong_header_order, no_header_meta):
    
            with self.assertRaises(ValueError):
                with bedio.open_param(ds.filename, "r") as fpar:
                    pass


    def test_no_header(self):

        no_header = sim_bed.ParamBedDataSet(os.path.join(bedio_dir.name,
                                            "no_header.bed"),
                                            header=None)

        with self.assertRaises(ValueError):
            with bedio.open_param(no_header.filename, "r") as fpar:
                pass

    def test_empty_file(self):
        empty_file = sim_bed.ParamBedDataSet(os.path.join(bedio_dir.name,
                                                          "empty_file.bed"),
                                             meta=None,
                                             header=None,
                                             data=None)

        with self.assertRaises(UnboundLocalError):
            with bedio.open_param(empty_file.filename, "r") as fpar:
                pass

    def test_not_ordered_gene_id(self):

        wrong_order = sim_bed.ParamBedDataSet(os.path.join(bedio_dir.name,
                                                           "unsorted_records.bed"),
                                              data = [sim_bed.DEFAULT_DATA[-1],
                                                      *sim_bed.DEFAULT_DATA])

        with bedio.open_param(wrong_order.filename, "r") as fpar:

            with self.assertRaises(ValueError):

                for key, vals in fpar.group_by("gene_id"):
                    pass

    def test_no_name(self):
        with self.assertRaises(TypeError):
            tmp = bedio.ParseParamBed()



if __name__ == "__main__":
    unittest.addModuleCleanup(tearDownModule)
    unittest.main()
    unittest.doModuleCleanup()
