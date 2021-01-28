# rcca

An R package for Sparse Canonical Correlation Analysis (CCA).

An implementation of the sparse CCA method proposed by Suo et al. (2017) extended to produce multiple canonical vector pairs (Rodosthenous et al 2020).

## Installation

### With `devtools`

```
library(devtools)
devtools::install_github("mkomod/rcca")
```

### From source

```
$ git clone https://github.com/mkomod/rcca
$ R CMD build rcca
$ R CMD INSTALL rcca_0.0.0.9000.tar.gz
```

*Note:* the version number `_0.0.0.9000` might be different

## References

Suo, X., Minden, V., Nelson, B., Tibshirani, R., & Saunders, M. (2017). Sparse canonical correlation analysis. Machine Learning, 83(3), 331–353. https://doi.org/10.1007/s10994-010-5222-7

Rodosthenous, T., Shahrezaei, V., & Evangelou, M. (2020). Integrating multi-OMICS data through sparse canonical correlation analysis for the prediction of complex traits: a comparison study. Bioinformatics, 36(17), 4616–4625. https://doi.org/10.1093/bioinformatics/btaa530
