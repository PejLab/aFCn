"""Create ficticious genotype data.

Genotype data is stored as VCF file, and used for 
makeing simulation data.

By: GDML
"""

from collections import OrderedDict
import re

import numpy as np


class VCFMetaInfoABC:
    _line_prefix = "##"
    _key_val_delim = "="

class VCFStructureMetaInfoLine(VCFMetaInfoABC):
    """VCF structured meta information line data.

    Args:
        key: (str) line key
        args: (str) with format 
            description attribute key = description attribute value,
            for example, 'ID=GT'
    """
    def __init__(self, key, *args):
        self._line_prefix = "##"
        self._key_val_delim = "="
        self._description_prefix = "<"
        self._description_suffix = ">"
        self._description_attribute_delim = ","

        self.key = key
        self.attributes = OrderedDict()

        for warg in args:
            key, val = warg.split(self._key_val_delim)

            self.attributes[key] = val

    def to_string(self):
        """Print vcf column descriptor meta information line."""
        return (self._line_prefix + self.key + 
                self._key_val_delim + self.description())

    def description(self):
        """Make VCF attribute description string.

        For example: <ID=test,Description=This is a simple example>
        """
        attributes = []

        for key, val in self.attributes.items():
            attributes.append(key + self._key_val_delim + val)

        return (self._description_prefix +
                self._description_attribute_delim.join(attributes) +
                self._description_suffix)


class VCFMetaInfoLine(VCFMetaInfoABC):
    def __init__(self, key, val):
        self.key = key
        self.val = val

    def to_string(self):
        return (self._line_prefix +
                self.key +
                self._key_val_delim +
                self.val)
    

class VCFRecord:
    _col_delim = "\t"
    _phase_delim = "|"
    _geno_entry_regex = "^[0-9.][|/][0-9.]:?\\w*$"
    _geno_phased_regex = "^[0-9.]\\|[0-9.]:?\\w*$"
    _geno_biallelic_regex = "^[01.][|/][01.]:?\\w*$"
    _geno_entry_cap_group = "^([0-9.])[|/]([0-9.]):?\\w*$"

    def __init__(self, chrom, pos, id_, ref, alt, 
                 qual, filt, info, format_,
                 samples):
        self.chrom = chrom
        self.pos = pos
        self.id = id_
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filt
        self.info = info
        self.format = format_
        self.samples = samples

    @property
    def genotypes(self):
        genos = np.zeros(shape=(2, len(self.samples)))

        for i, samp in enumerate(self.samples):
            tmp = re.findall(self._geno_entry_cap_group, samp)[0]

            genos[0, i] = int(tmp[0]) if tmp[0] != "." else np.nan
            genos[1, i] = int(tmp[1]) if tmp[1] != "." else np.nan

        return genos

    def is_phased(self):
        for samp in self.samples:
            if re.match(self._geno_phased_regex, samp) is None:
                return False
        return True

#    def is_biallelic(self):
#        for samp in self.samples:
#            if re.match(self._geno_biallelic_regex, samp) is None:
#                return False
#
#        return True

    def is_qc_passed(self):
        return self.filter.lower() == "pass"

    def to_string(self):
        return self._col_delim.join([self.chrom,
                                     str(self.pos),
                                     self.id,
                                     self.ref,
                                     self.alt,
                                     str(self.qual),
                                     self.filter,
                                     self.info,
                                     self.format,
                                     *self.samples])


class VCFDataSet:
    _col_delim = "\t"
    _new_line = "\n"
    _header_prefix = "#"

    meta = [VCFMetaInfoLine("fileformat", "VCFv4.4"),
            VCFMetaInfoLine("fileDate","20230622"),
            VCFMetaInfoLine("source","No Source, this is a 'dummy' VCF"),
            VCFMetaInfoLine("reference","No reference"),
            VCFMetaInfoLine("phasing","Mixed")]

    structured_meta = [VCFStructureMetaInfoLine("contig",
                                                "ID=chrm10",
                                                "length=12345678",
                                                "assembly=no_assembly",
                                                "species=not_a_species"),
                       VCFStructureMetaInfoLine("INFO",
                                                "ID=SD",
                                                "Number=1",
                                                "Type=Integer",
                                                "Description=\"Ficticious Data\""),
                       VCFStructureMetaInfoLine("INFO",
                                                "ID=RD",
                                                "Number=1",
                                                "Type=Integer",
                                                "Description=\"Space holder\""),
                       VCFStructureMetaInfoLine("FILTER",
                                                "ID=q10",
                                                "Description=\"Qualtiy score below 10\""),
                       VCFStructureMetaInfoLine("FILTER",
                                                "ID=s50",
                                                "Description=\"Less than 50% of samples with genotype.\""),
                       VCFStructureMetaInfoLine("FORMAT",
                                                "ID=GT",
                                                "Number=1",
                                                "Type=String",
                                                "Description=\"Genotype\""),
                       VCFStructureMetaInfoLine("FORMAT",
                                                "ID=GQ",
                                                "Number=1",
                                                "Type=Integer",
                                                "Description=\"Genotype Quality\"")]
    samples = ["SAMPLE_1", "SAMPLE_2", "SAMPLE_3"]

    header = ["#CHROM","POS","ID",
              "REF","ALT","QUAL",
              "FILTER","INFO","FORMAT"] + samples

    records = [VCFRecord("chrm10",100,"var_1", "A", "T", 30, "PASS", 
                         "SD=10;RD=42", "GT:GQ", 
                        ["0|0:34", "0|1:20","1|0:21"]),
              VCFRecord("chrm10",1010,"var_2", "A", "T,C", 30, "PASS", 
                        "SD=12;RD=52", "GT:GQ", 
                        ["0|0:34", "0|2:20","1|2:21"]),
              VCFRecord("chrm10",101000,"var_3", "A", ".", 3, "q10", 
                        "SD=12;RD=52", "GT:GQ", 
                        ["0|0:34", "0|1:20","1|0:21"]),
              VCFRecord("chrm10",101002,"var_4", "T", "A", 21, "PASS", 
                        "SD=12;RD=52", "GT:GQ", 
                        ["0|0:34", "1|1:20","1|0:21"]),
              VCFRecord("chrm10",101003,"var_5", "T", "A", 21, "PASS", 
                        "SD=12;RD=52", "GT:GQ", 
                        ["0|0:34", "1/1:20","1|0:21"]),
              VCFRecord("chrm10",101003,"var_5", "T", "A", 21, "PASS", 
                        "SD=12;RD=52", "GT:GQ", 
                        ["0|0:34", "1/1:20","1|0:21"]),
              VCFRecord("chrm10",101004,"var_none_0", "T", "A", 21, "PASS", 
                        "SD=12;RD=52", "GT:GQ", 
                        ["0|0:34", "1/1:20",".|0:21"]),
              VCFRecord("chrm10",101007,"var_none", "A", "T", 21, "PASS", 
                        "SD=12;RD=52", "GT:GQ", 
                        ["0|0:34", ".|.:20","1|0:21"]),
              VCFRecord("chrm10",101008,"var_none_2", "A", "T", 21, "PASS", 
                        "SD=12;RD=52", "GT:GQ", 
                        ["0|0:34", "0|.:20","1|0:21"]),
              VCFRecord("chrm10",101012,"var_6", "T", "A", 1, "q10;s50", 
                        "SD=12;RD=52", "GT:GQ", 
                        ["0|0:34", "1|1:20","1|0:21"]),
              VCFRecord("chrm10",101103,"var_7", "T", "A", 21, "q10", 
                        "SD=12;RD=52", "GT:GQ", 
                        ["0|0:34", "1/1:20","1|0:21"]),
              VCFRecord("chrm10",101105,"var_8", "T", "A", 21, ".", 
                        "SD=12;RD=52", "GT:GQ", 
                        ["0|0:34", "1|1:20","1|0:21"])]

    def __str__(self):
        output_string = ""
        for meta_info in self.meta:
            output_string += meta_info.to_string() + self._new_line

        for struct_meta_info in self.structured_meta:
            output_string += (struct_meta_info.to_string() +
                              self._new_line)

        output_string += self._col_delim.join(self.header) + self._new_line

        for record in self.records:
            output_string += record.to_string() + self._new_line

        return output_string.strip()
