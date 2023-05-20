# Inferring allele Fold Change (aFC) from phased data

Mohammadi et al. 2017 [1] defined allele fold change (aFC) as the 
ratio in the number of gene transcripts under the alternative 
allele with respect to that of the reference allele.  Consequently,
it is a parameter that quantifies the effect of any one regulatory
variant with its target gene.  The authors mathematical showed how
to combine individual phased genotypes and aFC to predict observed
gene expression.  While this definition and model is general, the
authors used it to specifically study *cis*-regulatory effects 
and RNA sequencing data.  

Let's revisit the model by defining variables:

* $Y\in\mathbb{R_{>=0}}$ the abundance of a gene, 
* $x_h\in\{0,1\}^{n\times 1}$ encoding $n$ biallelic 
*cis*-regulatory, reference (0) and alternative (1) allele, variants for haplotype $h\in\{1,2\}$, 
* $\alpha\in\mathbb{R}$ the log abundance under the reference haplotype
* $\beta\in\mathbb{R}^{n\times 1}$ the log allele fold change

then it follows that

$$
\log(Y) = 
    \log\left(e^{\alpha + x_1^{T}\beta} + e^{\alpha + x_2^{T}\beta}\right) + 
    \epsilon
$$ 

where $\epsilon\sim\mathcal{N}(0,\sigma^2)$.

This software package for the Python programming language can be 
used to:

* **fit** log afc ($\beta$) from gene expression data,
* **predict** gene expression abundances from genotype data, and

from either the command line or within a Python script.

## Installation

The package can be installed directly from the GitHub repo using `pip`

```
python -m pip install git+
```

or cloned and installed from local source code using pip.



## Examples


## API


## References

[1] Mohammadi et al. *Genome Research* 2017

```
@article{Mohammadi2017GenomeResearch,
  title={Quantifying the regulatory effect size of cis-acting genetic variation using allelic fold change},
  author={Mohammadi, Pejman and Castel, Stephane E and Brown, Andrew A and Lappalainen, Tuuli},
  journal={Genome research},
  volume={27},
  number={11},
  pages={1872--1884},
  year={2017},
  publisher={Cold Spring Harbor Lab}
}
```

[2] Ehsan et al. BioRxiv

```
@article {Ehsan2022BioRxiv,
	author = {Nava Ehsan and Bence M. Kotis and Stephane E. Castel and Eric J. Song and Nicholas Mancuso and Pejman Mohammadi},
	title = {Haplotype-aware modeling of cis-regulatory effects highlights the gaps remaining in eQTL data},
	year = {2022},
	doi = {10.1101/2022.01.28.478116},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2022/01/28/2022.01.28.478116},
	journal = {bioRxiv}
}
```
