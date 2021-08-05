# This script should be executed in an R session after creating a conda environment
library("devtools")
download.file("https://github.com/campbio/celda/archive/refs/tags/v1.7.3.zip", destfile = "v1.7.3.zip")
devtools::install_local("v1.7.3.zip") # when prompted to update packages - update none
# for vision (including anndata for converting h5ad files to R objects):
install_version("anndata", version = "0.7.5.2", repos = "http://cran.us.r-project.org")
library("anndata")
anndata::install_anndata()
install_github("YosefLab/VISION")