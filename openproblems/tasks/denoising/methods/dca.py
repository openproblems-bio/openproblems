from ....tools.decorators import method
from ....tools.utils import check_version


def _dca(adata, test=False, epochs=None):
    from dca.api import dca

    import anndata

    if test:
        epochs = epochs or 30
    else:  # pragma: nocover
        epochs = epochs or 300

    # make adata object with train counts
    adata_train = anndata.AnnData(adata.obsm["train"])
    # run DCA
    dca(adata_train, epochs=epochs)

    # set denoised to Xmat
    adata.obsm["denoised"] = adata_train.X
    # check version of dca
    adata.uns["method_code_version"] = check_version("dca")
    return adata


@method(
    method_name="DCA",
    method_summary=(
        "DCA (Deep Count Autoencoder) is a method to remove the effect of dropout in"
        " scRNA-seq data. DCA takes into account the count structure, overdispersed"
        " nature and sparsity of scRNA-seq datatypes using a deep autoencoder with a"
        " zero-inflated negative binomial (ZINB) loss. The autoencoder is then applied"
        " to the dataset, where the mean of the fitted negative binomial distributions"
        " is used to fill each entry of the imputed matrix."
    ),
    paper_name="Single-cell RNA-seq denoising using a deep count autoencoder",
    paper_reference="eraslan2019single",
    paper_year=2019,
    code_url="https://github.com/theislab/dca",
    image="openproblems-python-tensorflow",
)
def dca(adata, test=False, epochs=None):
    return _dca(adata, test=test, epochs=epochs)
