from collections.abc import Iterable
import numpy as np


def sample_maf(j_snps: int,
               maf_min: float, maf_max: float,
               rng: np.random.Generator) -> np.ndarray[float]:
    """Sample minor allele frequence (maf)."""
    if (maf_max < maf_min 
        or maf_max >= 1 or maf_min <= 0):
        raise ValueError("maf values must be in (0,1) and maf_max > maf_min")

    if (j_snps <= 0):
        raise ValueError("Number of snps must be a positive integer")

    delta = maf_max - maf_min

    return maf_min + delta * rng.random(size=j_snps)


def sample_haplotypes(n_samples: int, maf: Iterable[float],
                      rng: np.random.Generator) -> tuple[np.ndarray]:
    """Sample biallelic haplotypes given maf. """

    j_snps = len(maf)
    hap_one = np.zeros(shape=(n_samples, j_snps))
    hap_two = np.zeros(shape=(n_samples, j_snps))

    for j, maf_val in enumerate(maf):
        hap_one[rng.random(size=n_samples) <= maf_val, j] = 1
        hap_two[rng.random(size=n_samples) <= maf_val, j] = 1

    return (hap_one, hap_two)