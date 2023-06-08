# HierMultinom
An R package implementing the method proposed in [Multiresolution categorical response regression for interpretable cell type annotation](https://arxiv.org/abs/2208.13857) by Molstad, A. J., and Motwani, K.

# Installation
MCMVR can be loaded directly into R through the the `devtools` package:
```{r}
install.packages("devtools")
library(devtools)
devtools::install_github("ajmolstad/HierMultinom")
```
# Citation
Please contact `amolstad@umn.edu` for citation directions. 

# Usage directions
Please visit [this example page](https://ajmolstad.github.io/docs/HierMultinomExample.html) for details on implementation and usage. This example is available above as an RMD file (HierMultinomExample.Rmd). 

# Simulation study results
One can reproduce the simulation study results using the scripts in the Simulations directory. Please see the README file in this directory for specifics. 

# scRNA-seq data analysis results and figures
One can reproduce the real single-cell RNA-seq data analyses using the scripts in the Application directory. Please see the README file in this directory for specifics. 
