# About
The package contains scripts aiming for dynamic tolerance peak matching and false discovery estimation for metabolite identification using spectral library search. The novelties of the package are:

1. Use peak matching instead of peak binning to calculate the cosine dot product score.
2. Dynamic mass tolerance in peak matching. We set dynamic tolerance on *m/z* for peak matching proportional with peak width according to the mass spectrometry type, which is proportional to (m/z)^2  for FTICR, to (m/z)^1.5 for Orbitrap, to *m/z* for Q-TOF, and is a constant for quadrupole mass analyzer.
3. Calculated the decoy score as the average of the cosine similarity score between the target and mass-shifted query spectra. 
4. Enable the calculation of Xcorr score based on the difference between target score and decoy score.
5. Enable the PEP score and FDR estimation based on the target and decoy strategy.

For more information about the dynamic theory, see the [article](https://www.sciencedirect.com/science/article/pii/S0003267021005006?via%3Dihub)

# System requirement

The package was tested with R (version 4.1.2) on a laptop equipped with an Intel Core i7-8550U CPU, 1 TB HDD, and 16GB RAM

# Installation

You can download the development version from [GitHub](https://github.com/xiaodfeng/DynamicTol) with:

```
if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("xiaodfeng/DynamicTol")
```
In case the installation of DynamicTol package is failed, you can load all the packages and functions at once by using the scripts in the Alternative way to use the package.R file of the vignettes folder.

# Bug reporting

Suggestions and bug reports are more than welcome at:<https://github.com/xiaodfeng/DynamicTol/issues>

## Citation

Please cite this package as: (To be filled in)

