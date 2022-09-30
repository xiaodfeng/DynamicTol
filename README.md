# About
The package contains scripts aiming for dynamic tolerance peak matching and false discovery estimation for metabolite identification using spectral library search. The novelties of the package are:

1. Use peak matching instead of peak binning to calculate the cosine dot product score.
2. Dynamic mass tolerance in peak matching. We set dynamic tolerance on *m/z* for peak matching proportional with peak width according to the mass spectrometry type, which is proportional to ![img](file:///C:/Users/xiaodong/AppData/Local/Temp/msohtmlclip1/01/clip_image002.gif) for FTICR, to ![img](file:///C:/Users/xiaodong/AppData/Local/Temp/msohtmlclip1/01/clip_image004.gif) for Orbitrap, to *m/z* for Q-TOF, and is a constant for quadrupole mass analyzer.
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

# Bug reporting

Suggestions and bug reports are more than welcome at:<https://github.com/xiaodfeng/DynamicTol/issues>

## Citation

Please cite this package as: (To be filled in)

