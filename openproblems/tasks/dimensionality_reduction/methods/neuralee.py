from ....tools.decorators import method
from ....tools.utils import check_version
from anndata import AnnData


@method(
    method_name="NeuralEE (CPU)",
    paper_name=" NeuralEE: A GPU-Accelerated Elastic Embedding "
    "Dimensionality Reduction Method for "
    "Visualizing Large-Scale scRNA-Seq Data ",
    paper_url="https://pubmed.ncbi.nlm.nih.gov/33193561/",
    paper_year=2020,
    code_url="https://github.com/HiBearME/NeuralEE",
    code_version=check_version("neuralee"),
    image="openproblems-python-extras",
)
def neuralee(adata: AnnData) -> AnnData:
    from neuralee.dataset import GeneExpressionDataset
    from neuralee.embedding import NeuralEE

    import torch

    dataset = GeneExpressionDataset(adata.X)
    # this follows "recommended" preprocessing taken from example notebooks, e.g.:
    # https://github.com/HiBearME/NeuralEE/blob/master/tests/notebooks/retina_dataset.ipynb
    dataset.log_shift()  # log(x + 1)
    dataset.subsample_genes(500)  # filter genes according to variance
    dataset.standardscale()
    # 10% of cells as batch to the estimate affinity matrix
    dataset.affinity_split(N_small=0.1)

    NEE = NeuralEE(dataset, d=2, device=torch.device("cpu"))
    res = NEE.fine_tune(verbose=False)

    adata.obsm["X_emb"] = res["X"]

    return adata
