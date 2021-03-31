from ....tools.decorators import method


def _marge(tss_to_peaks, adata):
    """https://genome.cshlp.org/content/26/10/1417.long"""
    import numpy as np
    import scipy

    alpha = -np.log(1.0 / 3.0) * 1e5 / 1e4
    tss_to_peaks["distance"] = tss_to_peaks.apply(
        lambda x: abs((int(x[5]) + int(x[6])) / 2 - int(x[1])) * 1.0 / 1e5, axis=1
    )
    decay_score = np.exp(-alpha * tss_to_peaks["distance"].values)
    tss_to_peaks["weight"] = 2.0 * decay_score / (1.0 + decay_score)
    gene_peak_weight = scipy.sparse.csr_matrix(
        (
            tss_to_peaks.weight.values,
            (tss_to_peaks.gene_index.astype("int32").values, tss_to_peaks.peak_index),
        ),
        shape=(adata.shape[1], adata.uns["mode2_var"].shape[0]),
    )
    adata.obsm["gene_score"] = adata.obsm["mode2"] @ gene_peak_weight.T


@method(
    method_name="MARGE",
    paper_name="""Modeling cis-regulation with a compendium of
                   genome-wide histone H3K27ac profiles""",
    paper_url="https://genome.cshlp.org/content/26/10/1417.long",
    paper_year=2016,
    code_version="1.0",
    code_url="https://github.com/suwangbio/MARGE",
    image="openproblems-python-extras",
)
def marge(adata, n_top_genes=500):
    from .beta import _atac_genes_score

    adata = _atac_genes_score(adata, top_genes=n_top_genes, method="marge")
    return adata
