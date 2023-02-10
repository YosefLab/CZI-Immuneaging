# This script should be executed in an R session after creating a conda environment
library("devtools")
download.file("https://github.com/campbio/celda/archive/refs/tags/v1.7.3.zip", destfile = "v1.7.3.zip")
devtools::install_local("v1.7.3.zip") # when prompted to update packages - update none