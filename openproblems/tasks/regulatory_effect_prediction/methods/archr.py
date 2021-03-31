from ....tools.decorators import method


def _archr_model21(tss_to_peaks, adata):
    """https://www.archrproject.com/bookdown/calculating-gene-scores-in-archr.html"""
    import numpy as np
    import scipy

    tss_to_peaks["distance"] = tss_to_peaks.apply(
        lambda x: abs((int(x[5]) + int(x[6])) / 2 - int(x[1])) * 1.0 / 5000, axis=1
    )
    tss_to_peaks["weight"] = np.exp(tss_to_peaks["distance"].values)
    gene_peak_weight = scipy.sparse.csr_matrix(
        (
            tss_to_peaks.weight.values,
            (tss_to_peaks.gene_index.astype("int32").values, tss_to_peaks.peak_index),
        ),
        shape=(adata.shape[1], adata.uns["mode2_var"].shape[0]),
    )
    adata.obsm["gene_score"] = adata.obsm["mode2"] @ gene_peak_weight.T


def _archr_model42(tss_to_peaks, adata):
    return


@method(
    method_name="ArchR",
    paper_name="""ArchR: An integrative and scalable software
                   package for single-cell chromatin accessibility analysis""",
    paper_url="https://www.biorxiv.org/content/10.1101/2020.04.28.066498v1",
    paper_year=2020,
    code_version="1.0",
    code_url="https://github.com/GreenleafLab/ArchR",
    image="openproblems-python-extras",
)
def archr_model21(adata, n_top_genes=500):
    from .beta import _atac_genes_score

    adata = _atac_genes_score(adata, top_genes=n_top_genes, method="archr_model21")
    return adata
