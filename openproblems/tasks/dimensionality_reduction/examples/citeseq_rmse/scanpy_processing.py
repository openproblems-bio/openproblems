import scanpy as sc


def run_standard_scanpy_preprocessing(adata):

    sc.pl.highest_expr_genes(
        adata,
        n_top=20,
    )
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    adata.var["mt"] = adata.var_names.str.startswith(
        "MT-"
    )  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
    )
    sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt")
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")

    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]

    sc.pp.normalize_total(adata, target_sum=1e4)

    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pl.highly_variable_genes(adata)

    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
    sc.pp.scale(adata, max_value=10)

    return adata


def scanpy_plot_pca(adata, color, edges, frameon, size):
    sc.pl.pca(adata, color=color, edges=edges, frameon=frameon, size=size)


def scanpy_plot_umap(adata, color, edges, frameon, size):
    sc.pl.umap(adata, color=color, edges=edges, frameon=frameon, size=size)


def scanpy_plot_tsne(adata, color, edges, frameon, size):
    sc.pl.tsne(adata, color=color, edges=edges, frameon=frameon, size=size)


def scanpy_pca(adata, color, edges, frameon, size):
    try:
        scanpy_plot_pca(adata, color, edges, frameon, size)
    except KeyError:
        sc.tl.pca(adata)
        sc.pl.pca_variance_ratio(adata, log=True)
        scanpy_plot_pca(adata, color, edges, frameon, size)


def scanpy_umap(adata, color, edges, frameon, size):
    try:
        scanpy_plot_umap(adata, color, edges, frameon, size)
    except KeyError:
        sc.tl.umap(adata)
        scanpy_plot_umap(adata, color, edges, frameon, size)


def scanpy_tsne(adata, color, edges, frameon, size):
    try:
        scanpy_plot_tsne(adata, color, edges, frameon, size)
    except KeyError:
        sc.tl.tsne(adata)
        scanpy_plot_tsne(adata, color, edges, frameon, size)


def scanpy_projection(
    adata,
    color=["leiden", "HUMAN_HDAC1", "HUMAN_TP53"],
    pca=True,
    umap=True,
    tsne=True,
    edges=False,
    frameon=False,
    size=30,
):
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
    sc.tl.leiden(adata)

    if pca:
        scanpy_pca(adata, color, edges, frameon, size)

    if tsne:
        scanpy_tsne(adata, color, edges, frameon, size)

    if umap:
        scanpy_umap(adata, color, edges, frameon, size)
