# Scale Simulation
Examples for "Scale Reliant Inference from Compositional Survey"

This repository contains all the code necessary to replicate the figures in *Scale Reliant Inference from Compositional Surveys* by Nixon, Letourneau, David, Mukherjee, and Silverman (2021+).

To run, first make sure the correct packages are installed. While most of packages are available on CRAN and/or Bioconductor, the `driver` package is available only on Github. To install:

```{r}
library(devtools)
devtools::install_github("jsilve24/driver")
```

To rerun the results, first make sure that you have the packages in the first few lines of  `02_simulations.R` and `03_realAnalysis.R` installed. If you are unsure of how to install packages from BioConductor (such as the \texttt{DESeq2} and \texttt{ALDEx2} package), see the instructions [here](https://www.bioconductor.org/install/).

Once the packages are installed, replication is easy! Simply run the functions in `01a_main_functions.R` and `01b_helper_functions.R`. Both the simulations and real data analysis files can then be run. They do not depend on each other, so run order does not matter.
