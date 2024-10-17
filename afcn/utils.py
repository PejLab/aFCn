"""Utility functions.

By: Genomic Data Modeling Lab
"""

# import datetime
# import os
# import numbers
from importlib import metadata
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
    """Test whether data set is biallelic (0,1, np.nan) numeric."""
    data_set = np.unique(x, equal_nan=True)

    if data_set.size > 3 or not is_numeric_nparray(data_set):
        return False

    delta = np.setdiff1d(data_set, np.array([0,1]))

    return delta.size == 0 or np.isnan(delta).all()


def get_version(file=None):
    return metadata.version("afcn")


def is_int(string_value):
    """Check whether input string is an integer."""
    return re.match("^[+-]?\\d+$", string_value) is not None


def is_float(string_value):
    """Check whether input string is a floating point number."""
    return re.match("^[+-]?\\d+\\.\\d*$", string_value) is not None


class Gradient:

    _h0 = 1e-3

    def __init__(self,f, tol=0.00000001):
        self._f = f
        self._tol = tol

    def __call__(self, x):
        p = x.size
        grad_vector = np.zeros(p)
        tmp = x.copy()

        for i in range(p):
            h = self._h0
            delta = 100 * self._tol
            dfdx = 1.

            while np.abs(delta > self._tol):
                previous_dfdx = dfdx

                tmp[i] = x[i] + h
                dfdx = self._f(tmp)

                tmp[i] = x[i] - h
                dfdx -= self._f(tmp)

                dfdx /= (2*h)

                delta = dfdx - previous_dfdx
                h /= 10

            tmp[i] = x[i]
            grad_vector[i] = dfdx

        return grad_vector


class HessianObject:
    def __init__(self, p):
        self.matrix = np.zeros(shape=(p, p))
        self.converged = None
        self.h = np.full(int(p * (p+1) / 2), -1)
        self.delta = np.full(int(p * (p+1) / 2), -1)



class Hessian:
    """Numerical computation of the second derivative.

    The calculation is adpative, changing the step size for
    which differences are computed.  Step size is selected
    such that the change in derivative between successive
    steps is less than `tol` or the minimum step size is 
    reached.

    Args:
        f (callable, e.g. function)
            A callable whose argument is only the position
            in which the Hessian is computed.
        tol: (int)
            Convergence criterion, the maximum difference between 
            successive difference calculations.
        
            
    """
    _h0 = 1e-3

    def __init__(self, f, tol = 1e-6, 
                 h_min=1e-12, h_fc=2):
        self._f = f
        self._tol = tol
        self._h_min = h_min
        self._h_fc = h_fc

    def __call__(self, x):
        p = x.size
        out = HessianObject(p)
    
        # compute diagonal entries
        tmp = x.copy()
        idx = 0
    
        for i in range(p):
            
            h = self._h0 * self._h_fc
            delta = self._tol * 100
            hess_ii = 1.
    
            while (np.abs(delta) > self._tol
                   and h >= self._h_min):

                h /= self._h_fc

                previous_hess_ii = hess_ii
                tmp[i] = x[i] + 2*h
                hess_ii = self._f(tmp)
        
                tmp[i] = x[i] - 2*h
                hess_ii += self._f(tmp)
    
                tmp[i] = x[i]
                hess_ii -= 2*self._f(tmp)
        
                hess_ii /= (4*(h**2))
        
                delta = hess_ii - previous_hess_ii
    
            out.matrix[i,i] = hess_ii
            out.h[idx] = h
            out.delta[idx] = np.abs(delta)

            if h < self._h_min:
                out.converged = False
                return out

            idx += 1
    
    
        tmp = x.copy()
        # compute the off-diagonal entries
        for i in range(p):
            for j in range(p):
    
                if i == j:
                    continue
    
                if i > j:
                    out.matrix[i,j] = out.matrix[j,i]
                    continue
    
                h = self._h0 * self._h_fc
                delta = self._tol * 100
                hess_ij = 1.
    
                while (np.abs(delta) > self._tol
                       and h >= self._h_min):

                    h /= self._h_fc

                    previous_hess_ij = hess_ij
    
                    tmp[i] = x[i] + h
                    tmp[j] = x[j] + h
                    hess_ij = self._f(tmp)
    
                    tmp[i] = x[i] + h
                    tmp[j] = x[j] - h
                    hess_ij -= self._f(tmp)
    
                    tmp[i] = x[i] - h
                    tmp[j] = x[j] + h
                    hess_ij -= self._f(tmp)
    
                    tmp[i] = x[i] - h
                    tmp[j] = x[j] - h
                    hess_ij += self._f(tmp)
    
                    hess_ij /= (4*h**2)
    
                    delta = hess_ij - previous_hess_ij
    
    
                out.matrix[i, j] = hess_ij
                out.h[idx] = h
                out.delta[idx] = np.abs(delta)

                if h < self._h_min:
                    out.converged = False
                    return out

                tmp[i] = x[i]
                tmp[j] = x[j]
                idx += 1

        if out.converged is None:
            out.converged = True

        return out
