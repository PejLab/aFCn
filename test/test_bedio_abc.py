""""Test for bedio abstract base classes

By: Robert Vogel
    Genomic Data Modeling Lab
"""
import os
from datetime import datetime
import tempfile
import io
from collections import OrderedDict
from unittest import TestCase, main

from afcn import bedio, utils



def setUpModule():

    global bed_dir
    global parse_bed_fname
    global write_bed_fname


    bed_dir = tempfile.TemporaryDirectory()
    parse_bed_fname = os.path.join(bed_dir.name, "parsing_test.bed")
    write_bed_fname = os.path.join(bed_dir.name, "write_test.bed")

    b = bedio.BedABC()

    s = ""
    for key, val in b.meta.items():
        s += (f"{b._meta_prefix}{key}{b._meta_data_key_val_delimiter}"
              f"{val}{b._new_line}")

    s += "".join([b._header_prefix,
                b._field_delimiter.join(["chrom","start", "end", "name"]),
                b._new_line])

    s += b._field_delimiter.join(["1", "1000", "1100", "ENS_TEST_001"])

    with open(parse_bed_fname, "w") as fid:
        fid.write(s)


def tearDownModule():
    bed_dir.cleanup()


class TestBedABC(TestCase):
    def setUp(self):
        _file_suffix = ".bed"
        _encoding = "utf-8"

        _meta_prefix = "##"
        _meta_data_key_val_delimiter = "="

        _header_prefix = "#"

        _field_delimiter = "\t"
        _new_line = "\n"

        _date = datetime.today()
        _date = f"{_date.year}-{_date.month:02d}-{_date.day:02d}"

        _py_version=".".join([str(os.sys.version_info.major),
                              str(os.sys.version_info.minor),
                              str(os.sys.version_info.micro)])

        if "USER" in os.environ:
            _user = os.environ["USER"]
        elif "LOGIN" in os.environ:
            _user = os.environ["LOGIN"]
        else:
            _user = "not_available"

        self.meta = OrderedDict(afcn_version=utils.get_version(),
                                date = _date,
                                python_version = _py_version,
                                user = _user)

    def test_bed_spec(self):
        b = bedio.BedABC()
        self.assertEqual(b._file_suffix , ".bed")

        self.assertEqual(b._encoding , "utf-8")

        self.assertEqual(b._meta_prefix , "##")
        self.assertEqual(b._meta_data_key_val_delimiter , "=")

        self.assertEqual(b._header_prefix , "#")

        self.assertEqual(b._field_delimiter , "\t")
        self.assertEqual(b._new_line , "\n")


    def test_initialize(self):
        b = bedio.BedABC()

        self.assertIsNone(b._fid)
        self.assertIsNone(b.header)

        meta_data_keys = ("afcn_version", "date", "python_version", "user")

        for key in meta_data_keys:
            self.assertIn(key, self.meta)

        self.assertEqual(len(self.meta), len(meta_data_keys))

        self.assertIsInstance(b.meta, OrderedDict)




class TestParseBedABC(TestCase):
    def setUp(self):
        self._fid = io.FileIO(parse_bed_fname, "w")
        self.parser = bedio.ParseBedABC(self._fid)

    def tearDown(self):
        if not self._fid.closed:
            self._fid.close()

    def test_initialize(self):
        self.assertEqual(self._fid, self.parser._fid)

        self.assertIsInstance(self.parser._colname_to_idx, dict)
        self._req_header_fields = None

    def test_context_management(self):
        self.assertFalse(self.parser.closed)
        self.parser.close()
        self.assertTrue(self.parser.closed)
        self.assertTrue(self._fid.closed)

        with bedio.ParseBedABC(io.FileIO(parse_bed_fname, 'r')) as fid:
            self.assertFalse(fid.closed)

        self.assertTrue(fid.closed)



if __name__ == "__main__":
    unittest.addModuleCleanup(tearDownModule)
    main()
    unittest.doModuleCleanup()
