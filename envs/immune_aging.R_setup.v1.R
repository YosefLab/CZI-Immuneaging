# run this file using Rscript
library("devtools")
install_version("anndata", version = "0.7.5.2", repos = "http://cran.us.r-project.org")
library("anndata")
anndata::install_anndata()
install_github("YosefLab/VISION")