"""Utility functions.

By: Genomic Data Modeling Lab
"""

import numbers
import os
import re
import numpy as np

NUMPY_NUMERIC_DTYPE_KINDS = ("f", "u", "i")


def is_numeric_nparray(x):
    """Test whether numpy array is numeric.

    Note: that a mixture of booleans and numeric may be

    """
    return x.dtype.kind in NUMPY_NUMERIC_DTYPE_KINDS


# TODO what to do about np.nan
# TODO is it OK that object arrays raise TypeError as opposed to
# returning false, see unit test for example

def is_biallelic(x):
    """Test whether data set is biallelic (0,1) numeric."""
    data_set = np.unique(x, equal_nan=True)

    if data_set.size > 2 or not is_numeric_nparray(data_set):
        return False

    return np.setdiff1d(data_set, 
                        np.array([0,1])).size == 0


def get_version(file=None):

    if file is None:
        file = os.path.join(os.path.dirname(__file__), "__init__.py")

    found_version_line = False

    with open(file, "r") as fid:

        for tline in fid:
            if re.search("__version__", tline) is not None:
                found_version_line = True
                break
            
    if (not found_version_line or 
        (tmp := re.search("\w+\.\w+\.\w+", tline)) is None):
        return None

    return tmp.group()

