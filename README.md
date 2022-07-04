# SpectralUnmixR

Unmix spectral flow cytometry data in R.

The following unmixing algorithms are currently implemented:

- Ordinary least squares (OLS)
- Weighted least squares (WLS)
- Non-negative least squares (NNLS)
- Weighted non-negative least squares (WNNLS)
- Mean absolute percentage error (MAPE)

Helper functions for the calculation of the reference spectra and interactive plot functions for the in depth-visualization of the unmixing results are included.


## Installation

First install flowCore from Bioconductor, then SpectralUnmixR and its dependencies:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("flowCore"))

if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("stbenke/SpectralUnmixR")
```

## Licence
This R package is licensed under the GPL-3 license (see LICENSE file).
