""""Tests for the prediction bed file IO.

By: GDML
"""
from unittest import TestCase, main
from datetime import datetime
import io

from afcn import bedio, utils

class TestWritePredictionBed(TestCase):
    def setUp(self):
        self.encoding = "utf-8"

        # meta data
        date = datetime.today()
        meta_data = [f"##afcn_version={utils.get_version()}\n",
                    f"##date={date.year}-{date.month:02d}-{date.day:02d}\n"]

        self.b_id = io.BytesIO()
        for line in meta_data:
            self.b_id.write(line.encode(encoding=self.encoding))

        # header
        self.sample_names = ["sample_1", "sample_2"]
        header = "\t".join(["#chrom",
                            "start",
                            "end",
                            "gene_id",
                            *self.sample_names])
        header += "\n"
        self.b_id.write(header.encode(encoding=self.encoding))

        self.b_id.seek(0)
        # data records


    def tearDown(self):
        if not self.b_id.closed:
            self.b_id.close()

    def test_init(self):
        fid = bedio.WritePredictionBed(self.b_id)

        self.assertEqual(fid._fid, self.b_id)

        self.assertFalse(fid._meta_and_header_written)

        self.b_id.close()
        self.assertTrue(fid.closed)

    def test_context_management(self):
        fid = bedio.WritePredictionBed(self.b_id)

        self.assertFalse(fid.closed)
        self.assertFalse(self.b_id.closed)

        self.assertEqual(fid, fid.__enter__())
        fid.__exit__()

        self.assertTrue(fid.closed)
        self.assertTrue(self.b_id.closed)
        
    def test_write_meta_data(self):
        with bedio.WritePredictionBed(io.BytesIO()) as fid:

            fid.write_meta_data(self.sample_names)
            fid._fid.seek(0)

            for correct_line in self.b_id:
                if not correct_line.startswith(b"#"):
                    break

                self.assertEqual(fid._fid.readline(), correct_line)

            self.assertTrue(fid._meta_and_header_written)
        
    def test_exeception(self):
        with (bedio.WritePredictionBed(io.BytesIO()) as fid,
              self.assertRaises(ValueError)):
            fid.write_line_record("1", 10, 100, "SIM_GENE_1", [3, 4.3])


if __name__ == "__main__":
    main()
