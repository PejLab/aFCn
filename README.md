![unit-tests](https://github.com/PejLab/aFCn/actions/workflows/unit-tests.yml/badge.svg?branch=new_interface)


# `afcn` A tool to fit and predict gene expression by a mechanistic model

The `afcn` program applies a mechanistic model of gene
expression regulation by *cis*-regulatory elements developed
by [Mohammadi et al. 2017 (1)](README.md#(1)) and 
[Ehsan et al. 2024 (2)](README.md#(2)).  Here, we provide two
submodules

* `afcn fit` to infer model parameters from data
* `afcn predict` to generate gene expression predictions from
    phased genotypes.


## Installation

The package requires `Python >= 3.9` and can be installed directly
from the GitHub repo using `pip`

```
```

Alternatively, clone this repository then install from
local source code using `pip`, we hope to support PyPI soon.


## Examples

Examples and ficticious data can be found in the [example directory](afcn/example/)
of this repository.  These examples demonstrate the enumerate options and
how to use each submodule.  (`afcn fit` is under construction).


## Method description

Mohammadi et al. 2017 (1) defined allele fold change (aFC) as the 
ratio in the number of gene transcripts under the alternative 
allele with respect to that of the reference allele.  Consequently,
it is a parameter that quantifies the effect of any one regulatory
variant with its target gene.  The authors of (1) showed mathematically
how to combine individual phased genotypes and aFC values to predict 
observed gene expression.  While this definition and model is general,
the authors used it to specifically study *cis*-regulatory effects 
of gene regulation.  

This software package for the Python programming language can be 
used as:

`afcn fit` 🚧 **under construction** 🚧

Infers model parameters ($\alpha$, $\beta$) values by least squares.  As an example
    consider the *cis*-regulation of an arbitrary gene.  Let
    $i = 1,2,\dots, N$ be an index identifying one of $N$ samples,
    $j = 1,2, \dots, J$ be an index identifying one of $J$ *cis*-regulatory
    loci, and $h=1,2$ be an index identifying one of 2 phased haplotypes.
    For each $i$, we have RNA Sequencing derived gene counts $y_i$ and
    $J$ length vectors of phased haplotypes $x_i^{(1)}$ and 
    $x_i^{(2)}$.  An allele at locus $j$ for haplotype $h$ sample $i$,
    $x_{ij}^{(h)}$, takes values 0 and 1 representing the presence of
    either reference or alternative allele, respectively.  Let's denote the
    model predicted expression of our arbitrary gene by haplotype
    $h$, defined by Mohammadi et al. 2017 (1),
  
$$
    g\big(x_i^{(h)}, \alpha, \beta\big) = 2^{\alpha + x_{i}^{(h) T}\beta}
$$

then the total expression of the gene in sample $i$ is

$$
f\big(x_i^{(1)}, x_i^{(2)},\alpha, \beta\big) = g\big(x_i^{(1)}, \alpha, \beta\big) \
    + g\big(x_i^{(2)}, \alpha, \beta\big)
$$
    
Where $\alpha$ is a scalar represents the log2 reference expression
and $\beta$ is a $J$ length vector the log2 fold change per locus.  Inference
of these parameters will be computed by least squares

$$
\hat{\alpha},\hat{\beta} = \underset{\alpha,\beta}{\text{argmin}} \
    \sum_{i=1}^N \left( \
    \log_2\big(y_i + 1\big) - \log_2\big( f(x_i^{(1)},x_i^{(2)}, \alpha,\beta)\big) \
    \right)^2
$$

`afcn predict`

Estimate the gene count of sample $i$ attributed to
  haplotype $h$ is $g(x_{i}^{(h)},\alpha=0,\beta)$, and
  the total gene count $f(x_{i}^{(1)},x_{i}^{(2)},\alpha=0,\beta)$.
  

## References

### (1) 

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

### (2)

```
@article{Ehsan2024NatureCommunications,
  title={Haplotype-aware modeling of cis-regulatory effects highlights the gaps remaining in eQTL data},
  author={Ehsan, Nava and Kotis, Bence M and Castel, Stephane E and Song, Eric J and Mancuso, Nicholas and Mohammadi, Pejman},
  journal={Nature Communications},
  volume={15},
  number={1},
  pages={522},
  year={2024},
  publisher={Nature Publishing Group UK London}
}
```
