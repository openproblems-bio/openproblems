from ....tools.decorators import method
from ....tools.normalize import log_cpm_hvg
from ....tools.utils import check_version
from anndata import AnnData
from typing import Optional

import logging

log = logging.getLogger("openproblems")


def _create_neuralee_dataset(
    adata: AnnData,
    subsample_genes: Optional[int] = 500,
    normalize: bool = True,
    batch_size: int = 1000,
):
    """Create a neuralee GeneExpressionDataset from an AnnData

    This follows "recommended" preprocessing taken from example notebooks, e.g.:
    https://github.com/HiBearME/NeuralEE/blob/master/tests/notebooks/retina_dataset.ipynb
    """
    from neuralee.dataset import GeneExpressionDataset

    dataset = GeneExpressionDataset(adata.X)
    if normalize:
        dataset.log_shift()  # log(x + 1)
    if subsample_genes is not None:
        dataset.subsample_genes(subsample_genes)  # filter genes according to variance
    if normalize:
        dataset.standardscale()
    # 1000 cells as a batch to estimate the affinity matrix
    dataset.affinity_split(N_small=min(batch_size, adata.n_obs))
    return dataset


@method(
    method_name="NeuralEE (CPU) (Default)",
    paper_name=" NeuralEE: A GPU-Accelerated Elastic Embedding "
    "Dimensionality Reduction Method for "
    "Visualizing Large-Scale scRNA-Seq Data ",
    paper_url="https://www.frontiersin.org/articles/10.3389/fgene.2020.00786/full",
    paper_year=2020,
    code_url="https://github.com/HiBearME/NeuralEE",
    code_version=check_version("neuralee"),
    image="openproblems-python-extras",
)
def neuralee_default(adata: AnnData, test: bool = False) -> AnnData:
    from neuralee.embedding import NeuralEE

    import torch

    # Store raw counts for use by metrics
    adata.layers["counts"] = adata.X.copy()

    # this can fail due to sparseness of data; if so, retry with more genes
    # note that this is a deviation from the true default behavior, which fails
    # see https://github.com/openproblems-bio/openproblems/issues/375
    subsample_genes = 500
    while True:
        try:
            dataset = _create_neuralee_dataset(adata, subsample_genes=subsample_genes)
        except ValueError:
            if subsample_genes < adata.n_vars:
                subsample_genes = min(adata.n_vars, int(subsample_genes * 1.2))
                log.warning(
                    "ValueError in neuralee_default. "
                    f"Increased subsample_genes to {subsample_genes}"
                )
            else:
                raise
        else:
            break

    NEE = NeuralEE(dataset, d=2, device=torch.device("cpu"))
    res = NEE.fine_tune(verbose=False)

    adata.obsm["X_emb"] = res["X"].detach().cpu().numpy()

    return adata


@method(
    method_name="NeuralEE (CPU) (logCPM, 1kHVG)",
    paper_name=" NeuralEE: A GPU-Accelerated Elastic Embedding "
    "Dimensionality Reduction Method for "
    "Visualizing Large-Scale scRNA-Seq Data ",
    paper_url="https://pubmed.ncbi.nlm.nih.gov/33193561/",
    paper_year=2020,
    code_url="https://github.com/HiBearME/NeuralEE",
    code_version=check_version("neuralee"),
    image="openproblems-python-extras",
)
def neuralee_logCPM_1kHVG(adata: AnnData, test: bool = False) -> AnnData:
    from neuralee.embedding import NeuralEE

    import torch

    adata = log_cpm_hvg(adata)

    dataset = _create_neuralee_dataset(adata, normalize=False, subsample_genes=None)

    NEE = NeuralEE(dataset, d=2, device=torch.device("cpu"))
    res = NEE.fine_tune(verbose=False)

    adata.obsm["X_emb"] = res["X"].detach().cpu().numpy()

    return adata
