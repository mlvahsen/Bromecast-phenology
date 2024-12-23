# Bromecast common garden phenology manuscript

Code and data for "Phenological sensitivity of *Bromus tectorum* genotypes depends on current and source environments" (Vahsen et al., in press at *Ecology*)

All code was run in R version 4.3.2 (2023-10-31). Necessary packages to run script are identified at the top of each script.

All R scripts are written assuming the below organizational structure:
```
Bromecast-phenology
│   README.md
└─── code
└─── data
└─── figs
└─── outputs
└─── supp_code
└─── supp_data
```
The main analysis file ```code/lmm.R``` will need to be run in order to generate outputs needed for subsequent analyses, graphs, and tables. Output from the models (*i.e.*, the model objects) will be stored in ```outputs/``` which are then accessed by subsequent scripts.

[![DOI](https://zenodo.org/badge/639471014.svg)](https://doi.org/10.5281/zenodo.14417844)
