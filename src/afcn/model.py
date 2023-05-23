"""Implement gene expression prediction model and fitting.

By: Genomic Data Model Lab


Available Functions:

    predict : given model parameters and phased cis-regulation genotypes
        predict the abundances of one genes transcripts.


    simulate : generate simulated gene expression from genotype
        and model parameters
"""
import numbers
import numpy as np
from . import utils


# TODO how to handle missing genotype data, e.g. nans, should I ignore?
# TODO how to handle missing effect sizes?
def _predict(hap_one, hap_two, alpha, beta):
    """Compute model, without input checks.

    Args:
        see predict
    """
    return (np.exp2(alpha + np.dot(hap_one, beta)) + 
            np.exp2(alpha + np.dot(hap_two, beta)))


def predict(hap_one, hap_two, alpha, beta):
    """Expected abundance of a gene's transcripts under model.

    Apply the haplotype aware model of gene expression published
    in Mohammadi et al. Genome Research 2017.

    Args:
        hap_one: ((n variants,) ndarray) or 
            ((N samples, n variants) ndarray) 
            biallelic genotypes where 0 denotes reference and 1 
            denotes the alternative alleles of first haplotype.
        hap_two: genotypes for haplotype 2, see hap_one
        alpha: (float) the log reference expression
        beta: ((n variants,) ndarray) log allele fold change,
            should always be a 1-d array

    Returns:
        either:
            ((N samples,) ndarray) with > 1 sample and arbirary number
                of variants
    """
    if (not utils.is_biallelic(hap_one) or
        not utils.is_biallelic(hap_two)):
        raise ValueError("Genotypes are not biallelic")

    if beta.ndim != 1 or not utils.is_numeric_nparray(beta):
        raise ValueError("beta must be 1-D np.ndarray.")

    if not isinstance(alpha, numbers.Number):
        raise ValueError("alpha must be a number, int or float")

    if hap_one.ndim == 1:
        hap_one = hap_one.reshape(1, hap_one.size)
    if hap_two.ndim == 1:
        hap_two = hap_two.reshape(1, hap_two.size)

    return _predict(hap_one, hap_two, alpha, beta)


def simulate(hap_one, hap_two, alpha, beta, sd, seed=None):
    """Simulate gene expression data.
    
    Args:
        hap_one: ((n variants,) ndarray) or
            ((N samples , n variants) ndarray)
            biallelic genotypes where 0 denotes reference and 1
            denotes the alternative alleles of first haplotype.
        hap_two: genotypes for haplotype 2, see hap_one
        alpha: (float) the log reference expression
        beta: ((n variants,) ndarray) log allele fold change,
            should always be a 1-d array
        sd: (float) standard deviation of log-normal noise
        seed: seed for random number generations with
            numpy.random.defaul_rng(seed), (default None)

    Returns:

    """
    rng = np.random.default_rng(seed=seed)

    size = 1
    if hap_one.ndim == 2:
        size = hap_one.shape[0]

    return (predict(hap_one, hap_two, alpha, beta) * 
            np.exp(rng.normal(0, sd, size=size)))

