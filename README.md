## Soil Security Laboratory (University of Sydney)

![](https://img.shields.io/badge/release-v0.1-blue.svg)


This is (will be) a compilation of function used for Digital Soil Mapping, and pedometrics in general.

So far, the functions included are:
  
  - Fuzzy k-means with extragrades (de Gruijter and McBratney (1988) ,Tranter *et al.*, (2010))

### Installation

This package can be installed usign  [devtools](https://github.com/hadley/devtools).

```R
> devtools::install_github("spadarian/soilsecuritylab")
```

For instructions on how to install devtools, please refer to their decumentation.

### References

- De Gruijter, J. and McBratney, A. 1988. A modified fuzzy k-means method for predictive classification. In: 1st Conference of the International Federation of Classification Societies. Ed. by H. Bock. Elsevier Science, Amsterdam: pp. 97–104.
- Tranter, G, Minasny, B., and McBratney, A. 2010. Estimating Pedotransfer Function Prediction Limits Using Fuzzy k-Means with Extragrades. Soil Sci. Soc. Am. J. 74 (6): 1967–1975.

### Todos

- Finish implementation of equal-area quadratic splines.
- Improve documentation.

### License

GPL-3
