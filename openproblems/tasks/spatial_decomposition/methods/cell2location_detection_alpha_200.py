from ....tools.decorators import method
from ....tools.utils import check_version
from ..utils import split_sc_and_sp

import numpy as np


@method(
    method_name="Cell2location_detection_alpha_200",
    paper_name="Cell2location maps fine-grained cell types in spatial transcriptomics",
    paper_url="https://www.nature.com/articles/s41587-021-01139-4",
    paper_year=2022,
    code_url="https://github.com/BayraktarLab/cell2location",
    image="openproblems-python-extras",
)
def stereoscope(adata, test=False, max_epochs_sc=None, max_epochs_sp=None):

    from cell2location.models import Cell2location
    from cell2location.models import RegressionModel

    if test:
        max_epochs_sp = max_epochs_sp or 10
        max_epochs_sc = max_epochs_sc or 10
        num_samples = 10
    else:  # pragma: nocover
        max_epochs_sc = max_epochs_sc or 250
        max_epochs_sp = max_epochs_sp or 30000
        num_samples = 1000

    adata_sc, adata = split_sc_and_sp(adata)

    if "tech" in adata_sc.obs.columns:
        adata_sc.obs["batch_key"] = adata_sc.obs["tech"].copy()
    else:
        adata_sc.obs["batch_key"] = "all"

    # REFERENCE SIGNATURE ESTIMATION FROM SCRNA
    # prepare anndata for the regression model
    RegressionModel.setup_anndata(
        adata=adata_sc,
        # 10X reaction / sample / batch
        batch_key="batch_key",
        # cell type, covariate used for constructing signatures
        labels_key="label",
    )
    sc_model = RegressionModel(adata_sc)
    sc_model.train(max_epochs=max_epochs_sc)
    # In this section, we export the estimated cell abundance
    # (summary of the posterior distribution).
    adata_sc = sc_model.export_posterior(
        adata_sc, sample_kwargs={"num_samples": num_samples, "batch_size": 2500}
    )
    # export estimated expression in each cluster
    if "means_per_cluster_mu_fg" in adata_sc.varm.keys():
        inf_aver = adata_sc.varm["means_per_cluster_mu_fg"][
            [
                f"means_per_cluster_mu_fg_{i}"
                for i in adata_sc.uns["mod"]["factor_names"]
            ]
        ].copy()
    else:
        inf_aver = adata_sc.var[
            [
                f"means_per_cluster_mu_fg_{i}"
                for i in adata_sc.uns["mod"]["factor_names"]
            ]
        ].copy()
    inf_aver.columns = adata_sc.uns["mod"]["factor_names"]

    # SPATIAL MAPPING
    # find shared genes and subset both anndata and reference signatures
    intersect = np.intersect1d(adata.var_names, inf_aver.index)
    adata = adata[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    # prepare anndata for cell2location model
    adata.obs["sample"] = "all"
    Cell2location.setup_anndata(adata=adata, batch_key="sample")
    # create and train the model
    st_model = Cell2location(
        adata,
        cell_state_df=inf_aver,
        # the expected average cell abundance: tissue-dependent
        # hyper-prior which can be estimated from paired histology:
        # here = average in the simulated dataset
        N_cells_per_location=20,
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection:
        detection_alpha=200,
    )

    st_model.train(
        max_epochs=max_epochs_sp,
        # train using full data (batch_size=None)
        batch_size=None,
        # use all data points in training because
        # we need to estimate cell abundance at all locations
        train_size=1,
    )
    # In this section, we export the estimated cell abundance
    # (summary of the posterior distribution).
    adata = st_model.export_posterior(
        adata,
        sample_kwargs={
            "num_samples": num_samples,
            "batch_size": st_model.adata.n_obs,
        },
    )

    adata.obsm["proportions_pred"] = adata.obsm["q05_cell_abundance_w_sf"]
    adata.uns["method_code_version"] = check_version("cell2location")
    return adata
