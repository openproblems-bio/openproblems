from scanpy import AnnData

from ....tools.decorators import method


@method(
    method_name="No denoising",
    paper_name="N/A",
    # paper_url="https://books.google.com/books?id=64JYAwAAQBAJ",
    # paper_year=2013,
    # code_url="https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html",
    code_version=check_version("scanpy"),
)
def no_denoising(adata: AnnData) -> None:
    """Do nothing."""
    adata.obsm["denoised"] = adata.obsm["train"]
