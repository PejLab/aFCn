import io
import itertools
import numpy as np
import gzip
from pysam import VariantFile

# TODO how to manage specification version, e.g. code may 
# change but spec does not.
# TODO where to write specification so that it is callable 
# to __main__.py and this code


def read_bed(filename):
    """Read either gzipped or normal text file.

    Args:
        filename: (str) either .bed or .bed.gz file

    Returns:
        returns instance of class that mimics / uses io.FileIO
        services
    """

    if filename.endswith(".bed.gz"):
        return ParseParamGzipBed(filename)

    elif filename.endswith(".bed"):
        return ParseParamBed(filename)

    raise ValueError("Not a bed file.")


class ParseParamBedABC:
    """Enforce parameter bed file specification."""
    _meta_prefix = "##"
    _meta_data_key_val_delimiter = "="

    _header_prefix = "#"
    _req_header_fields = ("chrom", "start", "end", "qtl_id", 
                          "gene_start", "gene_end", "gene_id",
                          "variant_pos", "variant_id", 
                          "ref", "alt", "log2_afc")
    _colname_to_idx  = dict()

    _field_delimiter = "\t"

    def __init__(self):

        if not hasattr(self, "name"):
            raise NotImplementedError

        self.meta = dict()
        self.header = None

        # Note: if the initialization fails, the context manager methods
        # will not be called, in this case, close file object if
        # open
        self._initialize()

    def __iter__(self):
        raise NotImplementedError

    def __next__(self):
        raise NotImplementedError

    def _initialize(self):
        """Load and validate meta data and header as defined in spec."""
        for idx, col_name in enumerate(self._req_header_fields):
            self._colname_to_idx[col_name] = idx

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
        try:
            par_line = par_line.removeprefix(self._header_prefix)
        except UnboundLocalError as err:
            err.add_note(("\nMost likely reason for failure is that "
                          f"file {self.name} is empty."))
            raise 

        if not par_line.startswith(self._req_header_fields[0]):
            raise ValueError("Header must be present")

        self.header = par_line.strip().split(sep=self._field_delimiter)

        for i, field_name in enumerate(self._req_header_fields):

            if field_name != self.header[i]:
                raise ValueError(f"Invalid file header for {self.name}")

    def _record_parser(self):
        if self.header is None:
            return None

        for record in self:
            yield record.strip().split(self._field_delimiter)

    def idx(self, col_name):
        return self._colname_to_idx[col_name]

    def group_by(self, gb_id):
        """Perform group by on gb_id, assumed sorted.

        Requires that gene id are sorted lexicographically.
        """

        # start with ! as it has lowest lexicographic order
        previous_record_gb_id = "!"

        for k, g in itertools.groupby(self._record_parser(), 
                                      key=lambda x: x[self.idx(gb_id)]):

            if k < previous_record_gb_id:
                raise ValueError(f"{gb_id} column is not sorted.")

            previous_record_gb_id = k

            yield k, g


class ParseParamBed(io.FileIO, ParseParamBedABC):
    def __init__(self, filename):

        io.FileIO.__init__(self, filename)

        try:
            ParseParamBedABC.__init__(self)
        except:
            if not self.closed:
                self.close()

            raise

    def __next__(self):
        return io.FileIO.__next__(self).decode()


class ParseParamGzipBed(gzip.GzipFile, ParseParamBedABC):
    def __init__(self, filename):

        gzip.GzipFile.__init__(self, filename)

        try:
            ParseParamBedABC.__init__(self)
        except:
            if not self.closed:
                self.close()
                raise

    def __next__(self):
        return gzip.GzipFile.__next__(self).decode()


#def open_vcf(vcf):
#    if not os.exists(f"{vcf}.tbi"):
#        raise ValueError("Need tabix index for vcf file.")
#
#    return VariantFile(vcf, index_filename=f"{vcf}.tbi")
#
#
#class GenotypeBuffer:
#    def __init__(self, n_max_vars, n_samples):
#        self._buffer = np.zeros(shape=(n_max_vars, n_samples))
#        self._rows = 0
#        self._cols = 0
#
#    def reset(self):
#        self._rows = 0
#        self._cols = 0
#
#    @property
#    def size(self):
#        return self._rows * self._cols
#
#    @property
#    def shape(self):
#        return self._rows, self._cols
#
#
#class GenotypeVcfParser:
#    def __init__(fvcf, half_interval=500000):
#        self.filename = fvcf.filename.decode()
#        self._fid = fvcf
#        self.half_interval = cis_half_interval
#        self._n_samples = len(fvcf.header.samples)
#        self._hap_one_genotype = GenotypeBuffer(10000, self._nsamples)
#        self._hap_two_genotype = GenotypeBuffer(10000, self._nsamples)
#    
#    @property
#    def contig_len(self, name):
#        return self._fid.header.contigs[name].length
#
#    def get_variant_genotypes(self, chrom, pos):
#        # remember that fetch always returns a Tabix Iterator object
#        variant = next(self._fid.fetch(chrom, pos-1, pos))
#
#        
#        for i, samp_obj in enumerate(variant):
#            data[i]
#        return geneotype.allele_indices
#
#
#    def cis_variant_genotypes(self, chrom, tss):
#        
#        start = tss - self.half_interval
#
#        if start < 0:
#            start = 0
#
#        end = tss + self.half_interval
#
#        if end > self.contig_len[chrom]:
#            end = self.contig_len[chrom]
#
#        for variants in self._fid.fetch(chrom, start, end):
#            pass
#
#    


# class ParameterFileWriter
