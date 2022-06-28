from ....tools.decorators import method
from ....tools.utils import check_version
from .._utils import split_sc_and_sp


@method(
    method_name="Tangram",
    paper_name="Deep learning and alignment of spatially resolved single-cell transcriptomes with Tangram",  # noqa: E501
    paper_url="https://www.nature.com/articles/s41592-021-01264-7",
    paper_year=2021,
    code_url="https://github.com/broadinstitute/Tangram",
    image="openproblems-python-extras",
)
def tangram(adata, test=False, num_epochs=1000, n_markers=100):
    # analysis based on: https://github.com/broadinstitute/Tangram/blob/master/tutorial_tangram_with_squidpy.ipynb # noqa: E501
    # using tangram from PyPi, not github version

    # import dependencies
    import pandas as pd
    import scanpy as sc
    import tangram as tg
    import torch as t

    # get single cell reference
    ad_sc, adata = split_sc_and_sp(adata)
    # pre-process single cell data
    sc.pp.normalize_total(ad_sc, 1e4)
    sc.pp.log1p(ad_sc)
    # identify marker genes
    sc.tl.rank_genes_groups(ad_sc, groupby="label", use_raw=False)

    # extract marker genes to data frame
    markers_df = pd.DataFrame(ad_sc.uns["rank_genes_groups"]["names"]).iloc[
        0:n_markers, :
    ]

    # get union of all marker genes
    markers = list(set(markers_df.melt().value.values))

    # match genes between single cell and spatial data
    tg.pp_adatas(ad_sc, adata, genes=markers)

    # get device
    device = t.device("cuda" if t.cuda.is_available() else "cpu")

    # map single cells to spatial locations
    ad_map = tg.map_cells_to_space(
        ad_sc,
        adata,
        device=device,
        num_epochs=num_epochs,
    )

    # transfer labels from mapped cells to spatial location
    tg.project_cell_annotations(adata_map=ad_map, adata_sp=adata, annotation="label")

    # normalize scores
    pred_props = adata.obsm["tangram_ct_pred"].copy()
    pred_props_sum = pred_props.values.sum(axis=1, keepdims=True)
    adata.obsm["proportions_pred"] = pred_props.iloc[:, :] / pred_props_sum

    # remove un-normalized predictions
    del adata.obsm["tangram_ct_pred"]

    adata.uns["method_code_version"] = check_version("tangram-sc")

    return adata
