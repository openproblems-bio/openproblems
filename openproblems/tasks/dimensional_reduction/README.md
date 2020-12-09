## Dimensional Reduction

Dimensional reduction is a key challenge in single-cell data representation [<a href="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02128-7">1</a>].

## API

**Methods** should assign dimensionally-reduced 2D embedding coordinates to `adata.obsm['X_emb']` as well as the method of choice, i.e., `adata.obsm['X_pca']` or `adata.obsm['X_tsne']`.

**Metrics** should calculate the quality or "goodness of fit" of a dimensional reduction **method**.

1. Raimundo, F., Vallot, C. & Vert, J. Tuning parameters of dimensionality reduction methods for single-cell RNA-seq analysis. Genome Biol 21, 212 (2020). https://doi.org/10.1186/s13059-020-02128-7
