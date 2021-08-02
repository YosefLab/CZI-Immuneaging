# conda activate immune_aging.py_env.v2
library("VISION")
library("anndata")

args = commandArgs(trailingOnly=TRUE)
h5ad_file = args[1] # should have .h5ad suffix
outfile = args[2]

data <- read_h5ad(h5ad_file)
meta <- data$obs
X <- t(as.matrix(data$X))
latent <- data$obsm$X_scVI
rownames(latent) <- colnames(X)
vis <- Vision(X, meta = meta, latentSpace=latent)
umap <- data$obsm$X_umap_scvi
rownames(umap) <- colnames(X)
vis <- addProjection(vis, "UMAP", umap)
vis <- analyze(vis)
saveRDS(vis, paste0(outfile,".rds"))