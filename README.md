# F_Family_cova
Codes to compute the _F_-family of covariance functions.

This repository have the files required to compute and reproduce the results of the paper regarding the _F_-Family paper (that we can not cite yet). 

The data used for the paper is stored on the folder `data` whilst the folders `c` and `r` contains the C and R code needed to run the main files, respectively.

In order to compile the C files, the user must install the libraries GSL (https://www.gnu.org/software/gsl/) and ARB (http://arblib.org). After this, the user must compile the C files `matriz_cov.c` and `mat_pred.c` using the flags  `-lgsl -larb -lflint` as follows:

```
R CMD SHLIB ./c/matriz_cov.c  -lgsl -larb -lflint
R CMD SHLIB ./c/mat_pred.c  -lgsl -larb -lflint
```

The file `data_analysis.R` contains all the codes regarding the analysis of the data and allows to reproduce the results. Some code can have a long waiting time, so the results obtained by the authors were stored into the file `Data_results.Rdat`, which can be loaded as an R image using the command `load("Data_results.Rdat")`.

The file `Investigating_anisotropy.R` computes and reproduces the variogram plot that is shown on the paper. Because it can be time consuming, the needed objects to reproduce the plots were stored into the image file `Data_results.Rdat`.

The file `CreateFigures.R` reproduces the figures and the animations included in the paper. The user can see the animated plots in the folder `Animations`.

## Authors

* **No name displayed yet** - *Initial work*
