"""Utility functions.

By: Genomic Data Modeling Lab
"""

import numbers
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
