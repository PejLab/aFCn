"""IO tools for afcn bed files.

By: Robert Vogel
    Genomic Data Modeling Lab


Available Functions

    open_param
        returns a file handle for parsing or writing parameter bed file
    open_predict 
        returns a file handle for parsing or writing prediction bed file
    read_eqtl_map
        returns a file handle for parsing eqtl_map file
    read_gene_expression
        returns a file handle for parsing gene expression data

"""

import os
import io
import itertools
from collections import OrderedDict
import gzip
from datetime import datetime

import numpy as np

from . import utils

# TODO how to manage specification version, e.g. code may 
# change but spec does not.
# TODO where to write specification so that it is callable 
# to __main__.py and this code




class BedABC:
    """Define parameter and methods for I/O with all bed files.

    This class is not meant to be used by itself, but to be subclassed.  It
    contains properties defining: delimiters, encoding, file suffix, and
    prefixes for the three types of data (meta, header, records) contained in
    such documents.  Moreover, it includes methods for context managers (to
    open and close file handles), iteration, and finding traversing files.

    To subclass BedABC you'll need to call

        __init__

    and pass a file handle to

        _fid
        
    """
    _file_suffix = ".bed"
    _encoding = "utf-8"

    _meta_prefix = "##"
    _meta_data_key_val_delimiter = "="

    _header_prefix = "#"

    _hap_delimiter = "|"
    _field_delimiter = "\t"
    _new_line = "\n"

    def __init__(self):
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

        self._fid = None
        self.header = None

    def __iter__(self):
        return self

    def __next__(self):
        return self._fid.__next__().decode(self._encoding)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        if not self.closed:
            self.close()

    def tell(self):
        return self._fid.tell()

    def seek(self, *args):
        return self._fid.seek(*args)

    @property
    def closed(self):
        return self._fid.closed

    def close(self):
        if not self.closed:
            self._fid.flush()
            self._fid.close()


class WriteBedABC(BedABC):
    """Abstract base class for writing bed files.

    Use abstract base class by overloading:

        write_line_record (method)
            bed file specific method for writing as single record to file
        _req_header_fields (attribute)
            column names in header
    """
    def __init__(self, fid):

        super().__init__()

        self._fid = fid

        self._colname_to_idx  = dict()

        self._meta_and_header_written = False
        self._req_header_fields = None

    def _write_meta_and_header_data(self):

        if self._req_header_fields is None:
            raise AttributeError("Required header fields haven't been specfied.")

        if self._meta_and_header_written:
            raise ValueError("Meta and header data have already been written")

        # instantiate string
        s = ""

        # construct meta data
        for key, val in self.meta.items():
            s += "".join([self._meta_prefix,
                          key,
                          self._meta_data_key_val_delimiter,
                          val,
                          self._new_line])

        # construct header
        s += self._header_prefix
        s += self._field_delimiter.join(self._req_header_fields.keys())
        s = s.encode(encoding=self._encoding)

        if self._fid.write(s) != len(s):
            raise RuntimeError("Bytes written not equal to bytes of string.")

        self._meta_and_header_written = True

    def write_line_record(self, *args):
        raise NotImplementedError


class ParseBedABC(BedABC):
    """Abstract base class for parsing bed files
    
    Use abstract base class by overloading:

        _record_parser (method)
            bed file specific record parser
            
            
    """
    def __init__(self, fid):
        super().__init__()

        self._fid = fid

        self._colname_to_idx  = dict()
        self._req_header_fields = None

    def _initialize(self):
        """Load and validate meta data and header as defined in spec."""

        # load meta data 
        for par_line in self:

            if not par_line.startswith(self._meta_prefix):
                break

            # max_split kwarg in the string split method is important in 
            # the event that a user has an '=' in the value field.

            key, val = (par_line
                        .removeprefix(self._meta_prefix)
                        .split(sep=self._meta_data_key_val_delimiter,
                               maxsplit=1))

            self.meta[key.strip()] = val.strip()

        # decompose header and verify it is spec. compliant
        # remember that for empty files par_line is not associated 
        # with any value and throws an UnboundLocalError
        if not par_line.startswith(self._header_prefix):
            raise ValueError("Header required: No file header found")


        try:
            
            par_line = par_line.removeprefix(self._header_prefix)

        except UnboundLocalError as err:
            # the add_note method is available on Python 3.11 +
            #err.add_note(("\nMost likely reason for failure is that "
            #              f"file {self.name} is empty."))
            err.args = (err.args[0] + 
                        f"\nMost likely reason for failure is "
                        "that file is empty.",)
            raise UnboundLocalError(err) from None

        self.header = par_line.strip().split(sep=self._field_delimiter)

        for i, field_name in enumerate(self._req_header_fields.keys()):

            if field_name != self.header[i]:
                raise ValueError(f"Invalid file header")

        for i, colname in enumerate(self.header):
            self._colname_to_idx[colname] = i


        self._data_char_number = self.tell()

    def _record_parser(self):
        """Bed file specific record parser."""
        raise NotImplementedError

    def idx(self, col_name):
        """Retrieve index of a field specified header defined column name"""
        if col_name not in self._colname_to_idx:
            raise KeyError(f"{col_name} is not in header")
        return self._colname_to_idx[col_name]

    def group_by(self, gb_id):
        """Perform group by on gb_id, assumed sorted.

        Requires that gb_id are sorted lexicographically.
        """

        if (self.tell() != self._data_char_number):
            self.seek(self._data_char_number)

        # start with ! as it has lowest lexicographic order
        previous_record_gb_id = "!"

        for k, g in itertools.groupby(self._record_parser(), 
                                      key=lambda x: x[self.idx(gb_id)]):

            if k < previous_record_gb_id:
                raise ValueError(f"{gb_id} column is not sorted.")

            previous_record_gb_id = k

            yield k, list(g)


def generate_eqtl_req_fields():
    return OrderedDict(chrom = str,
                       qtl_start = int,
                       qtl_end = int,
                       qtl_id = str,
                       gene_start = int,
                       gene_end = int,
                       gene_id = str,
                       variant_pos = int,
                       variant_id = str,
                       ref = str,
                       alt = str)



def generate_param_bed_spec():
    """Enforce parameter bed file specification."""
    _req_header_fields = generate_eqtl_req_fields()

    _req_header_fields["log2_afc"] = float

    return _req_header_fields



def open_param(filename, mode):
    """Open either gzipped or normal text file for reading.

    Args:
        filename: (str)
            either .bed or .bed.gz file
        mode: (char)
            "r" for read, "w" for write

    Returns:
        returns instance of class that mimics / uses io.FileIO
        services
    """

    if (not filename.endswith(".bed.gz") and 
        not filename.endswith(".bed")):

        raise ValueError("Not a bed file.")

    if not os.path.exists(filename):
        raise FileNotFoundError(filename)

    if mode == "r" and filename.endswith(".bed"):
        return ParseParamBed(io.FileIO(filename, mode))

    elif mode == "r" and filename.endswith(".bed.gz"):
        return ParseParamBed(gzip.GzipFile(filename))

    if mode == "w" and filename.endswith(".bed"):
        return WriteParamBed(io.FileIO(filename, mode))

    raise NotImplementedError


class ParseParamBed(ParseBedABC):
    def __init__(self, fid):
        # recall that the order of class methods being called is:
        #   1. __init__
        #   2. __enter__
        #   3. __exit__
        # If there exists an error in __init__ or __enter__, the __exit__
        # method is not guaranteed to run.  Consequently, I use the try block.
        try:

            super().__init__(fid)
            self._req_header_fields = generate_param_bed_spec()
            self._initialize()

        except Exception:

            self.__exit__()
            raise

    def _record_parser(self):

        if self.header is None:
            return None

        for record in self:

            output = record.strip().split(self._field_delimiter)

            i = 0
            for _, field_type in self._req_header_fields.items():
                output[i] = field_type(output[i])
                i += 1

            for out_val in output[i:]:

                if utils.is_int(out_val):
                    out_val = int(out_val)
                elif utils.is_float(out_val):
                    out_val = float(out_val)

                output[i] = out_val

                i += 1

            yield output


# TODO
class WriteParamBed(WriteBedABC):
    def __init__(self, fid):
        super().__init__(fid)

        self._req_header_fields = generate_param_bed_spec()

        try:

            super().__init__(fid)

            self.meta["data"] = "Inferred afc model parameters."

        except Exception:

            self.__exit__()
            raise

    def _line_record(self):
        """Write line of predictions.

        Args:
            chrom: (str) chromosome name
            qtl_start: (int) genomic coordinate of gene beginning
            qtl_end: (int) genomic coordinate of gene end
            qtl_id: (str)
            gene_star: (int)
            gene_end: (int)
            gene_id: (int) Ensembl gene id
            variant_pos: (int) Single nucleotide polymorphism
            variant_id: (str),
            ref: (str), reference allele
            alt: (str), alternative allele
            log2_afc: (float)
            log2_afc_sem: (float)
            log2_afc_p_value: (float)


        Returns:
            None
        """
        pass
#         if not self._meta_and_header_written:
#             raise ValueError("Write meta data before real data")
# 
#         data_str = self._field_delimiter.join([str(w) for w in data])
#         chrom = self._new_line + chrom
# 
#         record_str = self._field_delimiter.join([chrom, 
#                                                 str(start), 
#                                                 str(end), 
#                                                 name,
#                                                 data_str])
#         record_str = record_str.encode(encoding=self._encoding)
# 
#         if self._fid.write(record_str) != len(record_str):
#             raise RuntimeError("Bytes written not equal to bytes of string.")


def read_eqtl_map(filename):
    if filename.endswith(".bed"):
        return ParseEqtlMap(io.FileIO(filename))

    if filename.endswith(".bed.gz"):
        return ParseEqtlMap(gzip.GzipFile(filename))

    raise NotImplementedError


# TODO

class ParseEqtlMap(ParseBedABC):
    def __init__(self, fid):
        super().__init__(fid)
        self._req_header_fields = generate_eqtl_req_fields()

        self._initialize()

    def _initialize(self):
        """Load meta data and header, create header to idx dictionary."""

        # load meta data 
        for par_line in self:

            if not par_line.startswith(self._meta_prefix):
                break

            # max_split kwarg in the string split method is important in 
            # the event that a user has an '=' in the value field.

            key, val = (par_line
                        .removeprefix(self._meta_prefix)
                        .split(sep=self._meta_data_key_val_delimiter,
                               maxsplit=1))

            self.meta[key.strip()] = val.strip()

        # decompose header and verify it is spec. compliant
        # remember that for empty files par_line is not associated 
        # with any value and throws an UnboundLocalError
        if not par_line.startswith(self._header_prefix):
            raise ValueError("Header required: No file header found")

        try:
            
            par_line = par_line.removeprefix(self._header_prefix)

        except UnboundLocalError as err:
            # the add_note method is available on Python 3.11 +
            # err.add_note(("\nMost likely reason for failure is that "
            #              f"file {self.name} is empty."))
            err.args = (err.args[0] +
                        f"\nNo header found, most likely reason for failure is"
                        "that file is empty.",)
            raise UnboundLocalError(err) from None

        self.header = par_line.strip().split(sep=self._field_delimiter)

        # hcol : header column name, rcol: required header column name
        for hcol, rcol in zip(self.header, self._req_header_fields.keys()):

            if hcol != rcol:
                raise ValueError(f"Invalid file header")

        self._data_char_number = self.tell()

    def group_by(self, ):
        pass
    

def read_gene_expression(filename):
    if filename.endswith(".bed"):
        return ParseGeneExpression(io.FileIO(filename))

    if filename.endswith(".bed.gz"):
        return ParseGeneExpression(gzip.GzipFile(filename))

    raise NotImplementedError


class ParseGeneExpression(ParseBedABC):
    def __init__(self, filename):
        raise NotImplementedError



def open_predict(filename, mode):
    """Open prediction file handle

    Args:
        filename (str)
        mode (char)
            read "r", or write "w"
    Returns:
        (WritePredictionBed) 
    """

    if mode == "r":
        raise NotImplementedError
    elif mode == "w":
        return WritePredictionBed(io.FileIO(filename, mode))

    raise ValueError("Uninterpretable mode")


def generate_prediction_bed_spec():
    """Prediction bed doc string."""
    return OrderedDict(chrom = str, 
                       start = int,
                       end = int,
                       gene_id = str)


class WritePredictionBed(WriteBedABC):
    """Buffered write of data to prediction bed file.

    Args:
        fid: (file object)
    """
    def __init__(self, fid):
        super().__init__(fid)

        self.meta["data"]="Predicted gene expression"
        
        self._req_header_fields = generate_prediction_bed_spec()

    def set_sample_names(self, sample_names):

        for samp_name in sample_names:
            self._req_header_fields[samp_name] = float

        self._write_meta_and_header_data()

    def write_line_record(self, chrom, start, end, name,
                          predictions, hap_two_predictions=None):
        """Write line of predictions.

        Args:
            chrom: (str) chromosome name
            start: (int) genomic coordinate of gene beginning
            end: (int) genomic coordinate of gene end
            name: (str) gene id
            predictions: ((n sample,) np.ndarray)
                the behavior changes according to the value of
                hap_two_predictions as follows:
                    = None then predictions are log2 (total gene expresssion)
                    = (n sample,) np.ndarray) then floats representing
                        predicted gene expression from haplotype one
            hap_two_predictions: ((n sample,) np.ndarray) or None
                of floats representing predicted gene expression
                from haplotype two.  (default is None)

        Return:
            None
        """
        if not self._meta_and_header_written:
            raise ValueError("Write meta data before real data")

        data_str = []

        if hap_two_predictions is None:

            for pexpr in predictions:

                data_str.append(str(pexpr))

        else:

            for h1, h2 in zip(predictions, hap_two_predictions):

                data_str.append(self._hap_delimiter.join([str(h1), str(h2)]))


        data_str = self._field_delimiter.join(data_str)
        chrom = self._new_line + chrom

        record_str = self._field_delimiter.join([chrom, 
                                                str(start), 
                                                str(end), 
                                                name,
                                                data_str])

        record_str = record_str.encode(encoding=self._encoding)

        if self._fid.write(record_str) != len(record_str):
            raise RuntimeError("Bytes written not equal to bytes of string.")
