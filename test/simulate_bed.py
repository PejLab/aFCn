import os
import gzip


TEST_DIR = os.path.join(os.path.dirname(__file__),
                        "tmp_test_data")


DEFAULT_META = dict(
                version = "1.21.30randomVersion",
                data_set = "random_numbers_for_testing_code"
                )

DEFAULT_HEADER = ("chrom","qtl_start","qtl_end","qtl_id",
                   "gene_start", "gene_end", "gene_id",
                   "variant_pos", "variant_id","ref", "alt",
                   "log2_afc", "sem","p_val")

DEFAULT_DATA = [("chr1",5000,100000,"SIMULATE_FEATURE_1/VAR_1",
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


class ParamBedDataSet:
    def __init__(self, filename,
                 meta = DEFAULT_META, header = DEFAULT_HEADER,
                 data = DEFAULT_DATA,
                 compress = False):                            

        self.filename = filename

        self.meta = meta    
        self.header= header 
        self.data = data

        self.write_mode = "w"

        self.compress = compress

        if self.compress:
            self.write_mode = "wb"

        self._fid = None
        self._write()

        self.col_2_idx = {}
        self.features = []

        if (self.header is not None
            and self.data is not None):

    
            for i, col_name in enumerate(self.header):
                self.col_2_idx[col_name] = i
    
            self.features = self._features()

    def _write(self):
        with open(self.filename, self.write_mode) as self._fid:
            self._fid.write(self.string)

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
            s += "#"
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
