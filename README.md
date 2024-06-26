# HierMultinom R Package
An R package implementing the method proposed in [Multiresolution categorical response regression for interpretable cell type annotation](https://arxiv.org/abs/2208.13857) by Molstad, A. J., and Motwani, K. (2023). Please contact amolstad@umn.edu with any questions or comments. 

Note that the package will be frequently updated. To create the results from the article exactly, please refer to the version of our software from June 10th, 2023. We hope that subsequent updates will improve computing times. 

### Installation
HierMultinom can be loaded directly into R through the `devtools` package:
```{r}
install.packages("devtools")
library(devtools)
devtools::install_github("ajmolstad/HierMultinom")
```
### Citation instructions
Please cite the most recent version of the article mentioned above. As of June 2023, this was the following (in bibtex): 
```
@article{molstad2023multiresolution,
  title={Multiresolution categorical regression for interpretable cell-type annotation},
  author={Molstad, Aaron J and Motwani, Keshav},
  journal={Biometrics},
  volume={79},
  number={4},
  pages={3485--3496},
  year={2023},
  publisher={Wiley Online Library}
}
```
### Vignette
Please visit [this example page](https://ajmolstad.github.io/docs/HierMultinomExample.html) for details on implementation and usage. 

### Simulation study results
The simulation study results can be reproduced using the scripts in the Simulations directory. Please see the README file in this directory for specifics. 

### scRNA-seq data analysis results and figures
The single-cell RNA-seq data analyses can be performed using the scripts in the Application directory. Please see the README file in this directory for specifics. 
