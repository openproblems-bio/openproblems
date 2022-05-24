# from ....tools.normalize import log_cpm
from .....tools.decorators import method
from .....tools.utils import check_version


@method(
    method_name="Combat",
    paper_name="Adjusting batch effects in microarray expression data using\
                empirical Bayes methods",
    paper_url="https://academic.oup.com/biostatistics/article/8/1/118/252073",
    paper_year=2007,
    code_url="https://scanpy.readthedocs.io/en/stable/api/scanpy.pp.combat.html",
    code_version=check_version("scanpy"),
    image="openproblems-python-batch-integration",  # only if required
)
def combat_full_unscaled(adata, test=False):
    from scib.integration import runCombat
    from scib.preprocessing import reduce_data

    adata = runCombat(adata, "batch")
    reduce_data(adata, umap=False)
    adata.obsm['X_emb'] = adata.obsm['X_pca']
    # Complete the result in-place
    return adata


@method(
    method_name="Combat (hvg/unscaled)",
    paper_name="Adjusting batch effects in microarray expression data using\
                empirical Bayes methods",
    paper_url="https://academic.oup.com/biostatistics/article/8/1/118/252073",
    paper_year=2007,
    code_url="https://scanpy.readthedocs.io/en/stable/api/scanpy.pp.combat.html",
    code_version=check_version("scanpy"),
    image="openproblems-python-batch-integration",  # only if required
)
def combat_hvg_unscaled(adata, test=False):
    from ._utils import hvg_batch
    from scib.integration import runCombat
    from scib.preprocessing import reduce_data

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = runCombat(adata, "batch")
    reduce_data(adata, umap=False)
    adata.obsm['X_emb'] = adata.obsm['X_pca']
    return adata


@method(
    method_name="Combat (hvg/scaled)",
    paper_name="Adjusting batch effects in microarray expression data using\
                empirical Bayes methods",
    paper_url="https://academic.oup.com/biostatistics/article/8/1/118/252073",
    paper_year=2007,
    code_url="https://scanpy.readthedocs.io/en/stable/api/scanpy.pp.combat.html",
    code_version=check_version("scanpy"),
    image="openproblems-python-batch-integration",  # only if required
)
def combat_hvg_scaled(adata, test=False):
    from ._utils import hvg_batch
    from ._utils import scale_batch
    from scib.integration import runCombat
    from scib.preprocessing import reduce_data

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    adata = runCombat(adata, "batch")
    reduce_data(adata, umap=False)
    adata.obsm['X_emb'] = adata.obsm['X_pca']
    return adata


@method(
    method_name="Combat (full/scaled)",
    paper_name="Adjusting batch effects in microarray expression data using\
                empirical Bayes methods",
    paper_url="https://academic.oup.com/biostatistics/article/8/1/118/252073",
    paper_year=2007,
    code_url="https://scanpy.readthedocs.io/en/stable/api/scanpy.pp.combat.html",
    code_version=check_version("scanpy"),
    image="openproblems-python-batch-integration",  # only if required
)
def combat_full_scaled(adata, test=False):
    from ._utils import scale_batch
    from scib.integration import runCombat
    from scib.preprocessing import reduce_data

    adata = scale_batch(adata, "batch")
    adata = runCombat(adata, "batch")
    reduce_data(adata, umap=False)
    adata.obsm['X_emb'] = adata.obsm['X_pca']
    return adata
