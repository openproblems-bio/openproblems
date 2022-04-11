from ....tools.decorators import method
from ....tools.normalize import preprocess_logCPM_1kHVG
from ....tools.utils import check_version
from anndata import AnnData


@method(
    method_name="NeuralEE (CPU) (Default)",
    paper_name=" NeuralEE: A GPU-Accelerated Elastic Embedding "
    "Dimensionality Reduction Method for "
    "Visualizing Large-Scale scRNA-Seq Data ",
    paper_url="https://pubmed.ncbi.nlm.nih.gov/33193561/",
    paper_year=2020,
    code_url="https://github.com/HiBearME/NeuralEE",
    code_version=check_version("neuralee"),
    image="openproblems-python-extras",
)
def neuralee_default(adata: AnnData, test: bool = False) -> AnnData:
    from neuralee.dataset import GeneExpressionDataset
    from neuralee.embedding import NeuralEE

    import torch

    dataset = GeneExpressionDataset(adata.X)
    # this follows "recommended" preprocessing taken from example notebooks, e.g.:
    # https://github.com/HiBearME/NeuralEE/blob/master/tests/notebooks/retina_dataset.ipynb
    dataset.log_shift()  # log(x + 1)
    dataset.subsample_genes(500)  # filter genes according to variance
    dataset.standardscale()
    # 1000 cells as a batch to estimate the affinity matrix
    dataset.affinity_split(N_small=min(1000, adata.n_obs))

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
    from neuralee.dataset import GeneExpressionDataset
    from neuralee.embedding import NeuralEE

    import torch

    adata = preprocess_logCPM_1kHVG(adata)

    dataset = GeneExpressionDataset(adata.X)
    dataset.affinity_split(N_small=min(1000, adata.n_obs))

    NEE = NeuralEE(dataset, d=2, device=torch.device("cpu"))
    res = NEE.fine_tune(verbose=False)

    adata.obsm["X_emb"] = res["X"].detach().cpu().numpy()

    return adata
