# HierMultinom
An R package implementing the method proposed in [Multiresolution categorical response regression for interpretable cell type annotation](https://arxiv.org/abs/2208.13857) by Molstad, A. J., and Motwani, K. (2023). Please contact amolstad@umn.edu with any questions or comments. 

# Installation
HierMultinom can be loaded directly into R through the the `devtools` package:
```{r}
install.packages("devtools")
library(devtools)
devtools::install_github("ajmolstad/HierMultinom")
```
# Citation instructions
Please cite the most recent version of the article mentioned above. As of June, 2023, this was the following (in bibtex). 
```
@misc{molstad2022multiresolution,
      title={Multiresolution categorical regression for interpretable cell type annotation}, 
      author={Aaron J. Molstad and Keshav Motwani},
      year={2022},
      eprint={2208.13857},
      archivePrefix={arXiv},
      primaryClass={stat.ME}
}
```
# Usage directions
Please visit [this example page](https://ajmolstad.github.io/docs/HierMultinomExample.html) for details on implementation and usage. This example is available above as an RMD file (HierMultinomExample.Rmd). 

# Simulation study results
The simulation study results can be reproduced using the scripts in the Simulations directory. Please see the README file in this directory for specifics. 

# scRNA-seq data analysis results and figures
The single-cell RNA-seq data analyses can be performed using the scripts in the Application directory. Please see the README file in this directory for specifics. 
