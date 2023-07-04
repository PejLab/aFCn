import os
import shutil
import io

import gzip
import unittest

from afcn import bedio

# TODO need to update so that temp data makes use of
# tempfile module in Python standard library


TEST_DIR = os.path.join(os.path.dirname(__file__),
                        "tmp_test_data")

SIM_META = dict(
    version = "1.21.30randomVersion",
    data_set = "random_numbers_for_testing_code"
)

SIM_COLS = ("chrom","qtl_start","qtl_end","qtl_id",
            "gene_start", "gene_end", "gene_id",
            "variant_pos", "variant_id","ref", "alt",
            "log2_afc","sem","p_val")

SIM_DATA = [
    ("chr1",5000,100000,"SIMULATE_FEATURE_1/VAR_1",
     7500, 100000, "SIMULATE_FEATURE_1",
     5000,"VAR_1","A","T", 2.2, 0.1, 0.002),
    ("chr1",5005, 100000, "SIMULATE_FEATURE_1/VAR_2",
     7500, 100000, "SIMULATE_FEATURE_1",
     5005, "VAR_2", "G", "T", 0.2, 0.1, 0.2),
    ("chr1",5007, 100000, "SIMULATE_FEATURE_1/VAR_3",
     7500, 100000, "SIMULATE_FEATURE_1",
     5007, "VAR_3", "G", "C", 0.01, 0.1, 0.5),
    ("chr2",10250, 1000000, "SIMULATE_FEATURE_2/VAR_4",
     10375, 100000, "SIMULATE_FEATURE_2",
     10250, "VAR_4", "A", "G", 0.1, 0.01, 0.0005)
]


# Use for in memory test of abstract base clase
class MockParser(io.BytesIO, bedio.ParseParamBedABC):
    def __init__(self, filename, data_string):

        self.name = filename

        io.BytesIO.__init__(self, data_string.encode(encoding="utf-8"))

        try:
            bedio.ParseParamBedABC.__init__(self)
        except:
            if not self.closed:
                self.close()

            raise

    def __next__(self):
        return io.BytesIO.__next__(self).decode()


def mk_test_fname(basename):
    return os.path.join(TEST_DIR, basename)


class DataSet:
    def __init__(self, filename, meta=None, header=None, data=None, 
                 compress=False):

        self.filename = filename
        self.meta = meta
        self.header = header
        self.data = data

        self.col_2_idx = {}
        self.features = []

        if header is not None:

            for i, col_name in enumerate(self.header):
                self.col_2_idx[col_name] = i

            self.features = self._features()

        self.write_mode = "w"

        self.compress = compress

        if self.compress:
            self.write_mode = "wb"


    def _features(self):
        feature = [self.data[0][self.col_2_idx["gene_id"]]]
        
        for record in self.data[1:]:

            tmp_feature = record[self.col_2_idx["gene_id"]]

            if tmp_feature != feature[-1]:
                feature.append(tmp_feature)

        return feature

    @property
    def string(self):

        s = ""

        if self.meta is not None:
            for key, val in self.meta.items():
                s += f"##{key} = {val}\n"
    
        if self.header is not None:
            s += "\t".join(self.header)
            s += "\n"
    
        if self.data is not None:
            for record in self.data:
                s += "\t".join([str(r) for r in record])
                s += "\n"
    
        s = s.strip()

        if self.compress:
            return gzip.compress(s.encode(encoding="utf-8"))

        return s


class FileIODataSet(DataSet):
    def write(self):
        with open(self.filename, self.write_mode) as fid:
            fid.write(self.string)


def setUpModule():
    # Only these data are written to disk and used for testing
    global FILE_IO_DATA
    FILE_IO_DATA = []

    # These data are in memory and make use of the MockParser
    # for testing the bedio.ParseParamBedABC

    global INVALID_HEADER
    INVALID_HEADER = []
    global NOT_ORDERED
    global EMPTY_FILE

    if not os.path.exists(TEST_DIR):
        os.mkdir(TEST_DIR)

    # write valid data sets
    FILE_IO_DATA.append(FileIODataSet(mk_test_fname("parameters.bed"),
                                      meta=SIM_META,
                                      header=SIM_COLS,
                                      data=SIM_DATA))

    FILE_IO_DATA.append(FileIODataSet(mk_test_fname("parameters.bed.gz"),
                                      meta=SIM_META,
                                      header=SIM_COLS,
                                      data=SIM_DATA,
                                      compress=True))

    for val in FILE_IO_DATA:
        val.write()


    INVALID_HEADER = [DataSet(mk_test_fname("no_header.bed"),
                              meta=SIM_META, 
                              data=SIM_DATA),
                      DataSet(mk_test_fname("no_header_meta_entries.bed"),
                              data=SIM_DATA),
                      DataSet(mk_test_fname("wrong_entries.bed"),
                              header=[SIM_COLS[1], SIM_COLS[0], *SIM_COLS[2:]],
                              data=SIM_DATA)]

    EMPTY_FILE = DataSet(mk_test_fname("empty_file.bed"))


    NOT_ORDERED = DataSet(mk_test_fname("feature_not_ordered.bed"),
                          meta=SIM_META,
                          header=SIM_COLS,
                          data = [SIM_DATA[-1],
                                  SIM_DATA[0]])



def tearDownModule():

    FILE_IO_DATA.clear()

    if os.path.exists(TEST_DIR):
        shutil.rmtree(TEST_DIR)


class Testopen_param(unittest.TestCase):

    def test_correct_inputs(self):
        pass
        #        for ds in FILE_IO_DATA:
        #
        #            #            if ds.filename.endswith(".bed.gz"):
        #            #                io_class_ = bedio.ParseParamGzipBed
        #            #            elif ds.filename.endswith(".bed"):
        #            #                io_class_ = bedio.ParseParamBed
        #
        #
        #            with bedio.open_param(ds.filename, "r") as tmp:
        #                self.assertIsInstance(tmp, io_class_)
    
    def test_not_bed_exception(self):
        with self.assertRaises(ValueError):
            bedio.open_param("test.bed.bz2", "r")

        with self.assertRaises(ValueError):
            bedio.open_param("test.vcf.gz", "r")

        with self.assertRaises(ValueError):
            bedio.open_param("bed.txt", "r")

        with self.assertRaises(ValueError):
            bedio.open_param("test_bed.gz", "r")


class TestParameterFileParser(unittest.TestCase):

    def test_file_build_fileio_data(self):

        for ds in FILE_IO_DATA:

            with bedio.open_param(ds.filename, "r") as fpar:

                for key, meta_val in ds.meta.items():
                    self.assertEqual(meta_val, fpar.meta[key])

                for i, header_val in enumerate(ds.header):
                    self.assertEqual(header_val, fpar.header[i])

    def test_group_by_fileio_data(self):

        for ds in FILE_IO_DATA:

            with bedio.open_param(ds.filename, "r") as fpar:

                feature_idx = 0
                record_idx = 0

                for feature, feature_data in fpar.group_by("gene_id"):
                    self.assertEqual(feature, ds.features[feature_idx])
                    feature_idx += 1

                    for var_feature_data in feature_data:

                        for true_val, parsed_val in zip(ds.data[record_idx], var_feature_data):

                            self.assertEqual(true_val, parsed_val)

                        record_idx += 1

    def test_invalid_headers(self):

        for data_set in INVALID_HEADER:

            with self.assertRaises(ValueError):
                with MockParser(data_set.filename, data_set.string) as fpar:
                    pass

    def test_empty_file(self):
        with self.assertRaises(UnboundLocalError):
            with MockParser(EMPTY_FILE.filename, EMPTY_FILE.string) as fpar:
                pass

    def test_not_ordered_gene_id(self):
        with MockParser(NOT_ORDERED.filename, NOT_ORDERED.string) as fpar:

            with self.assertRaises(ValueError):
                for key, vals in fpar.group_by("gene_id"):
                    pass

    def test_no_name(self):
        with self.assertRaises(AttributeError):
            tmp = bedio.ParseParamBedABC()



if __name__ == "__main__":
    unittest.addModuleCleanup(tearDownModule)
    unittest.main()
    unittest.doModuleCleanup()
