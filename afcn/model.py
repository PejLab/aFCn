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
from scipy import optimize as opt

from . import utils



# TODO how to handle missing genotype data, e.g. nans, should I ignore?
# TODO how to handle missing effect sizes?
def _predict(haplotype, alpha, beta):
    """Compute model, without input checks.

    Args:
        see predict
    """
    return np.exp2(alpha + np.dot(haplotype, beta))


def predict(haplotype, alpha, beta):
    """Expected abundance of a gene's transcripts under model.

    Apply the haplotype aware model of gene expression published
    in Mohammadi et al. Genome Research 2017.

    Args:
        haplotype: ((n variants,) ndarray) or ((N samples, n variants) ndarray) 
            biallelic genotypes where 0 denotes reference and 1 
            denotes the alternative alleles of first haplotype.
        alpha: (float)
            the log2 reference expression
        beta: ((n variants,) ndarray)
            log2 allele fold change, should always be a 1-d array

    Returns:
        (float,) or (N samples,) ndarray
            The abunance, linear scale, of transcripts predicted to
            originate from the input haplotype
    """
    if not utils.is_biallelic(haplotype):
        raise ValueError("Genotypes are not biallelic")

    if beta.ndim != 1 or not utils.is_numeric_nparray(beta):
        raise ValueError("beta must be 1-D np.ndarray.")

    if not isinstance(alpha, numbers.Number):
        raise ValueError("alpha must be a number, int or float")

    if haplotype.ndim == 1:
        haplotype = haplotype.reshape(1, haplotype.size)

    return _predict(haplotype, alpha, beta)


def simulate(hap_one, hap_two, alpha, beta, sd, seed=None):
    """Simulate gene expression data.
    
    Args:
        hap_one: ((n variants,) ndarray) or ((N samples , n variants) ndarray)
            biallelic genotypes where 0 denotes reference and 1
            denotes the alternative alleles of first haplotype.
        hap_two:
            genotypes for haplotype 2, see hap_one
        alpha: (float)
            the log2 reference expression
        beta: ((n variants,) ndarray)
            log2 allele fold change, should always be a 1-d array
        sd: (float)
            standard deviation of log-normal noise
        seed:
            seed for random number generations with
            numpy.random.defaul_rng(seed), (default None)

    Returns:

    """
    rng = np.random.default_rng(seed=seed)

    size = 1
    if hap_one.ndim == 2:
        size = hap_one.shape[0]

    return (predict(hap_one, alpha, beta)
                * np.exp(rng.normal(0, sd, size=size))
            + predict(hap_two, alpha, beta)
                * np.exp(rng.normal(0, sd, size=size)))


def _linear_expansion_model(haplotype_one, haplotype_two, y):
    """Linear regression using the Pseudo-inverse

    The linearized model is the first order Taylor expansion
    about all parameters equal to zero.  Under this linear
    model we can estimate initial parameter values to be used
    as seeds in the nonlinear model.

    Args:
        haplotype_one ((n_samples, j_variants) np.ndarray)
            biallelic genotypes
        haplotype_two ((n_samples, j_variants) np.ndarray)
            biallelic genotypes
        y ((n_samples,) np.ndarray)
            gene counts per sample

    Returns:
        (dict)
            "pars": parameters values
            "rank" : rank of linearized matrix X
    """
    n_samples = y.size
    genotypes = haplotype_one + haplotype_two
    j_variants = genotypes.shape[1]

    if genotypes.shape[0] != n_samples:
        raise ValueError("Genotypes and isoform abundance dims"
                         " do not agree")

    # data are the linearized model based features of
    # the genotype data.  Initialize with zero order term
    n_pars = (j_variants + 1) 

    # subtract zeroth order term
    y_effective = np.log(y + 1) - np.log(3)

    # Linearized feature matrix
    X = np.zeros(shape=(n_samples, n_pars))

    # alpha pars
    mult_factor = 1/6

    idx = 0
    X[:, idx] = 2 * mult_factor
    idx += 1
    for j in range(j_variants):
        X[:, idx] += genotypes[:, j] * mult_factor
        idx += 1

    omega = X.T @ X
    pars = np.linalg.pinv(omega) @ X.T @ y_effective

    return {
        "pars":pars,
        "rank":np.linalg.matrix_rank(omega),
        }



def _obj(haplotype_one, haplotype_two, y, reg, reg_const):
    """Function closure for least squares optimization.

    Return 
        function
    """

    if reg is not None:
        reg = reg.lower()

    if reg is None and reg_const is None:

        _penalty_func = lambda k: 0

    elif reg == 'l1' and isinstance(reg_const, numbers.Number):

        _penalty_func = lambda k: reg_const * np.sum(np.abs(k))

    elif reg == "l2" and isinstance(reg_const, numbers.Number):

        _penalty_func = lambda k: reg_const * k @ k

    else:
        raise ValueError(f"{reg} and {reg_const} are invalid"
                        " regularization method and constant pair")
    

    def _g(k):
        return (np.sum((y 
            - np.log2(1 + _predict(haplotype_one,
                               k[0],
                               k[1:])
                     + _predict(haplotype_two,
                                k[0],
                                k[1:])))**2)
            + _penalty_func(k))

    return _g


def fit(haplotype_one, haplotype_two, y, reg=None, reg_const=None):
    """Fit isoform abundances, + pseudo count, to isoform model.

    Fits data using Newton-CG minimizer

    Args:
        haplotype_one ((n_samples, j_variants) np.ndarray)
            biallelic genotypes
        haplotype_two ((n_samples, j_variants) np.ndarray)
            biallelic genotypes
        y ((n_samples,) np.ndarray)
            counts of each gene

    Returns:
       out 
    """
    # validate input for common mistakes
    if (not utils.is_biallelic(haplotype_one) 
        or not utils.is_biallelic(haplotype_two)):

        raise ValueError("Input haplotypes are not biallelic, e.g. (0,1,np.nan)")

    if haplotype_one.shape != haplotype_two.shape:
        raise ValueError("haplotypes are not identical dimension")

    if haplotype_one.shape[0] != y.size:
        raise ValueError("Number of samples in haplotype don't match predictions")

    if (y < 0).any():
        raise ValueError("Gene counts can't be negative")


    # find initial parameters values
    lout = _linear_expansion_model(haplotype_one,
                                   haplotype_two,
                                   y)

    # nonlinear optimization
    objective = _obj(haplotype_one, haplotype_two,
                     np.log2(y+1),
                     reg, reg_const)

    jacobian = utils.Gradient(objective)
    hessian = utils.Hessian(objective)

    out = opt.minimize(objective, 
                       lout["pars"],
                       method="Newton-CG",
                       jac = jacobian,
                       hess = hessian,
                       options={"disp":False})


    out.matrix_rank = lout["rank"]

    hessian = utils.Hessian(objective)
    out.hessian = hessian(out.x)
    return out
