from .....tools.decorators import method
from .....tools.utils import check_version


@method(
    method_name="scanvi",
    paper_name="Probabilistic harmonization and annotation of single‐cell\
                transcriptomics data with deep generative models",
    paper_url="https://www.embopress.org/doi/full/10.15252/msb.20209620",
    paper_year=2021,
    code_url="https://github.com/YosefLab/scvi-tools",
    code_version=check_version("scvi"),
    image="openproblems-python-batch-integration",  # only if required
)
def scanvi_full_unscaled(adata, test=False):
    from scib.integration import runScanvi
    from scanpy.preprocessing import neighbors

    adata.obs.rename(
        columns={"labels": "lab"}, inplace=True
    )  # ugly fix for scvi conversion error
    adata = runScanvi(adata, "batch", "lab")
    neighbors(adata, use_rep="X_emb")
    adata.obs.rename(
        columns={"lab": "labels"}, inplace=True
    )  # ugly fix for scvi conversion error
    return adata


@method(
    method_name="Scanvi (hvg/unscaled)",
    paper_name="Probabilistic harmonization and annotation of single‐cell\
                transcriptomics data with deep generative models",
    paper_url="https://www.embopress.org/doi/full/10.15252/msb.20209620",
    paper_year=2021,
    code_url="https://github.com/YosefLab/scvi-tools",
    code_version=check_version("scvi"),
    image="openproblems-python-batch-integration",  # only if required
)
def scanvi_hvg_unscaled(adata, test=False):
    from ._utils import hvg_batch
    from scib.integration import runScanvi
    from scanpy.preprocessing import neighbors

    adata.obs.rename(
        columns={"labels": "lab"}, inplace=True
    )  # ugly fix for scvi conversion error
    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = runScanvi(adata, "batch", "lab")
    neighbors(adata, use_rep="X_emb")
    adata.obs.rename(
        columns={"lab": "labels"}, inplace=True
    )  # ugly fix for scvi conversion error
    return adata


@method(
    method_name="Scanvi (hvg/scaled)",
    paper_name="Probabilistic harmonization and annotation of single‐cell\
                transcriptomics data with deep generative models",
    paper_url="https://www.embopress.org/doi/full/10.15252/msb.20209620",
    paper_year=2021,
    code_url="https://github.com/YosefLab/scvi-tools",
    code_version=check_version("scvi"),
    image="openproblems-python-batch-integration",  # only if required
)
def scanvi_hvg_scaled(adata, test=False):
    from ._utils import hvg_batch
    from ._utils import scale_batch
    from scib.integration import runScanvi
    from scanpy.preprocessing import neighbors

    adata.obs.rename(
        columns={"labels": "lab"}, inplace=True
    )  # ugly fix for scvi conversion error
    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    adata = runScanvi(adata, "batch", "lab")
    neighbors(adata, use_rep="X_emb")
    adata.obs.rename(
        columns={"lab": "labels"}, inplace=True
    )  # ugly fix for scvi conversion error
    return adata


@method(
    method_name="Scanvi (full/scaled)",
    paper_name="Probabilistic harmonization and annotation of single‐cell\
                transcriptomics data with deep generative models",
    paper_url="https://www.embopress.org/doi/full/10.15252/msb.20209620",
    paper_year=2021,
    code_url="https://github.com/YosefLab/scvi-tools",
    code_version=check_version("scvi"),
    image="openproblems-python-batch-integration",  # only if required
)
def scanvi_full_scaled(adata, test=False):
    from ._utils import scale_batch
    from scib.integration import runScanvi
    from scanpy.preprocessing import neighbors

    adata.obs.rename(
        columns={"labels": "lab"}, inplace=True
    )  # ugly fix for scvi conversion error
    adata = scale_batch(adata, "batch")
    adata = runScanvi(adata, "batch", "lab")
    neighbors(adata, use_rep="X_emb")
    adata.obs.rename(
        columns={"lab": "labels"}, inplace=True
    )  # ugly fix for scvi conversion error
    return adata
