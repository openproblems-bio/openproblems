from ....tools.decorators import method


@method(
    method_name="Cross-Magic-Imputation",
    image="openproblems-python-extras",
    paper_name=None,
    paper_url=None,
    paper_year=None,
    code_url=None,
    code_version="1.0",
)
def cross_magic(adata, n_steps=3, tol=1e-2):
    """The tool uses the neighborhhod of cells inferred from RNA-Seq to
    stabalize the ATAC signal through a weighted average among neighbors.
    It calculates the MAGIC kernel using the palantir
    diffusion maps. The kernel is applied to impute the ATAC signal.

    Parameters
    ----------
    adata: An AnnData Object with the ATAC peack counts in adata.obsm['mode2']
    n_steps: The number of times the smoothening kernel is applied.
    tol: The tolerance below which imputed counts are regarded as zeros.
    """

    import numpy as np
    import palantir
    import pandas as pd
    import scanpy as sc

    sc.pp.pca(adata, n_comps=50)
    dm_res = palantir.utils.run_diffusion_maps(
        pd.DataFrame(adata.obsm["X_pca"], index=adata.obs_names)
    )
    T_steps = dm_res["T"] ** n_steps
    T_steps = T_steps.astype(np.float32)
    imputed_atac = np.dot(T_steps, adata.obsm["mode2_noisy"])
    imputed_atac.data[imputed_atac.data < tol] = 0
    imputed_atac.eliminate_zeros()
    adata.obsm["mode2_denoised"] = imputed_atac
    return adata
