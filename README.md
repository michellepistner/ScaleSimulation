# Scale Reliant Inference

This repository contains all the code necessary to replicate the figures in *Scale Reliant Inference* by Nixon, Letourneau, David, Lazar, Mukherjee, and Silverman (2021+).

To run, first make sure the correct packages are installed. All packages are available on CRAN and/or Bioconductor. CRAN packages needed: *fido*, *tidyverse*, *gghighlight*, *ggpattern*, *cowplot*, *directlabels*, *matrixNormal*, and *LaplacesDemon*. Bioconductor packages needed: *phyloseq*, *DESeq2*, and *ALDEx2*.

To rerun the results, first make sure that you have the packages in the first few lines of  `02_simulations.R` and `03_realAnalysis.R` installed. If you are unsure of how to install packages from BioConductor (such as the *DESeq2* and *ALDEx2* package), see the instructions [here](https://www.bioconductor.org/install/).

Once the packages are installed, replication should be easy. Simply run the functions in `01a_main_functions.R` and `01b_helper_functions.R`. Both the simulations and real data analysis files can then be run with these two files alone. They do not depend on each other, so run order does not matter.
