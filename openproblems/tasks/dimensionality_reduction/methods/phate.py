from ....tools.decorators import method
from ....tools.normalize import log_cpm_hvg
from ....tools.normalize import sqrt_cpm
from ....tools.utils import check_version


def _phate(adata, test=False, n_pca=None):
    from phate import PHATE

    if test:
        n_pca = n_pca or 10
    else:  # pragma: no cover
        n_pca = n_pca or 100

    phate_op = PHATE(n_pca=n_pca, verbose=False, n_jobs=-1)
    adata.obsm["X_emb"] = phate_op.fit_transform(adata.X)
    return adata


@method(
    method_name="PHATE (default pre-processing)",
    paper_name="Visualizing Transitions and Structure for Biological Data Exploration",
    paper_url="https://www.nature.com/articles/s41587-019-0336-3",
    paper_year=2019,
    code_url="https://github.com/KrishnaswamyLab/PHATE/",
    code_version=check_version("phate"),
    image="openproblems-python-extras",
)
def phate_default(adata, test: bool = False, n_pca=None):
    adata = sqrt_cpm(adata)
    return _phate(adata, test=test, n_pca=n_pca)


@method(
    method_name="PHATE (logCPM, 1kHVG)",
    paper_name="Visualizing Transitions and Structure for Biological Data Exploration",
    paper_url="https://www.nature.com/articles/s41587-019-0336-3",
    paper_year=2019,
    code_url="https://github.com/KrishnaswamyLab/PHATE/",
    code_version=check_version("phate"),
    image="openproblems-python-extras",
)
def phate_logCPM_1kHVG(adata, test: bool = False, n_pca=None):
    adata = log_cpm_hvg(adata)
    return _phate(adata, n_pca=n_pca)
