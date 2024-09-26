""""Tests for the bedio classes that write to file

By: Robert Vogel
    Genomic Data Modeling Lab

"""
import os
from datetime import datetime
import tempfile
import io
import numpy as np
from collections import OrderedDict
from unittest import TestCase, main

from afcn import bedio, utils



def setUpModule():

    global bed_dir
    global fname

    bed_dir = tempfile.TemporaryDirectory()
    fname = os.path.join(bed_dir.name, "write_test.bed")


def tearDownModule():
    bed_dir.cleanup()


class TestWriteBedABC(TestCase):

    def setUp(self):
        self._fid = io.FileIO(fname, "w")
        self.writer = bedio.WriteBedABC(self._fid)

    def tearDown(self):
        if not self._fid.closed:
            self._fid.close()

    def test_initialize(self):
        self.assertEqual(self._fid, self.writer._fid)

        self.assertFalse(self.writer._meta_and_header_written)
        self.assertIsNone(self.writer._req_header_fields)

    def test_context_management(self):
        self.assertFalse(self.writer.closed)
        self.writer.__exit__()
        self.assertTrue(self.writer.closed)

    def test_required_header_fields(self):
        # required header fields are set by subclass, consequently
        # should raise an attribute error
        with self.assertRaises(AttributeError):
            self.writer._write_meta_and_header_data()

    def test_write_line_record(self):
        with self.assertRaises(NotImplementedError):
            self.writer.write_line_record()

    def test_protect_meta_data(self):
        self.writer._meta_and_header_written = True
        self.writer._req_header_fields = {"a":"b"}

        with self.assertRaises(ValueError):
            self.writer._write_meta_and_header_data()



class TestWritePrediction(TestCase):
    def setUp(self):
        self.encoding = "utf-8"

        _py_version=".".join([str(os.sys.version_info.major),
                              str(os.sys.version_info.minor),
                              str(os.sys.version_info.micro)])

        if "USER" in os.environ:
            user = os.environ["USER"]
        elif "LOGIN" in os.environ:
            user = os.environ["LOGIN"]
        else:
            user = "dont_know"

        # meta data
        date = datetime.today()
        meta_data = [f"##afcn_version={utils.get_version()}\n",
                     f"##date={date.year}-{date.month:02d}-{date.day:02d}\n",
                     f"##python_version={_py_version}\n",
                     f"##user={user}\n",
                     "##data=Predicted gene expression\n"]
        
        self.samp_names = ["samp1", "samp2", "samp3"]
        header = '\t'.join(["#chrom", "start", "end", "gene_id",
                  *self.samp_names])

        self.b_id = io.BytesIO()
        for line in meta_data:
            self.b_id.write(line.encode(encoding=self.encoding))

        self.b_id.write(header.encode(encoding=self.encoding))

        coord = 0
        rng = np.random.default_rng()
        self.loci_records = []
        self.n_samples = len(self.samp_names)
        self.n_records = 4
        for i in range(self.n_records):
            coord = 10*i + 1

            tmp = OrderedDict(chrom = "chrom1",
                              start=coord,
                              end=coord+1,
                              name=f"GENE_{i+1}")

            for key, val in tmp.items():
                if key == "chrom":
                    s =f"\n{val}"
                    continue

                s += f"\t{val}"

            tmp["data"] = np.zeros(self.n_samples)

            for j in range(self.n_samples):
                val = rng.normal()
                s += f"\t{val}"
                tmp["data"][j] = val

            self.loci_records.append(tmp)

            self.b_id.write(s.encode(encoding=self.encoding))

        # go to beginning of buffer representing the file we are writing to
        self.b_id.seek(0)


    def tearDown(self):
        if not self.b_id.closed:
            self.b_id.close()

    def test_write_meta_data(self):
        with bedio.WritePredictionBed(io.BytesIO()) as fid:

            self.assertFalse(fid._meta_and_header_written)
            fid.set_sample_names(self.samp_names)

            with self.assertRaises(ValueError):
                fid._write_meta_and_header_data()

            fid.seek(0)

            # a single '#' will account for meta data and headers
            for correct_line in self.b_id:
                if not correct_line.startswith(b"#"):
                    break

                self.assertEqual(fid._fid.readline().strip(),
                                 correct_line.strip())

            self.assertTrue(fid._meta_and_header_written)
      
    def test_exeception_unset_samples(self):
        with (bedio.WritePredictionBed(io.BytesIO()) as fid,
              self.assertRaises(ValueError)):

            fid.write_line_record("1", 10, 100, "SIM_GENE_1",
                                  [3, 4.3, 42])


    def test_write_line_record(self):

        with bedio.WritePredictionBed(io.BytesIO()) as fid:

            fid.set_sample_names(self.samp_names)
            for lr in self.loci_records:
                fid.write_line_record(lr["chrom"],
                                      lr["start"],lr["end"],
                                      lr["name"], lr["data"])

            fid.seek(0)

            for correct_val, test_val in zip(self.b_id, fid):
                self.assertEqual(correct_val.decode(encoding=self.encoding),
                                 test_val)





if __name__ == "__main__":
    unittest.addModuleCleanup(tearDownModule)
    main()
    unittest.doModuleCleanup()
