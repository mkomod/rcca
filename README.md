
# rcca

An R package for Sparse Canonical Correlation Analysis (CCA).

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

## Experimental

There is an experimental branch with multi-threaded functions for cross validating the hyperparameters.
