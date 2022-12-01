from ....tools.decorators import method
from ....tools.utils import check_version
from ..utils import split_sc_and_sp
from typing import Optional

import functools
import numpy as np

_cell2location_method = functools.partial(
    method,
    paper_name="Cell2location maps fine-grained cell types in spatial transcriptomics",
    paper_url="https://doi.org/10.1038/s41587-021-01139-4",
    paper_year=2022,
    code_url="https://github.com/BayraktarLab/cell2location",
    image="openproblems-python-extras",
)


def _cell2location(
    adata,
    detection_alpha,
    n_cells_per_location=20,
    hard_coded_reference=True,
    amortised=False,
    num_samples=None,
    sc_batch_size=2500,
    st_batch_size=None,
    test=False,
    max_epochs_sc=None,
    max_epochs_st=None,
):

    from cell2location.cluster_averages.cluster_averages import compute_cluster_averages
    from cell2location.models import Cell2location
    from cell2location.models import RegressionModel
    from torch.nn import ELU

    if test:
        max_epochs_sc = max_epochs_sc or 2
        max_epochs_st = max_epochs_st or 2
        num_samples = num_samples or 2
    else:  # pragma: nocover
        max_epochs_sc = max_epochs_sc or 250
        max_epochs_st = max_epochs_st or 30000
        num_samples = num_samples or 1000

    adata_sc, adata = split_sc_and_sp(adata)

    if not hard_coded_reference:
        # if "tech" in adata_sc.obs.columns:
        #     adata_sc.obs["batch_key"] = adata_sc.obs["tech"].copy()
        # else:
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
        sc_model.train(max_epochs=max_epochs_sc, batch_size=sc_batch_size)
        # In this section, we export the estimated cell abundance
        # (summary of the posterior distribution).
        adata_sc = sc_model.export_posterior(
            adata_sc,
            sample_kwargs={"num_samples": num_samples, "batch_size": sc_batch_size},
        )
        # export estimated expression in each cluster
        try:
            means_per_cluster = adata_sc.varm["means_per_cluster_mu_fg"]
        except KeyError:
            # sometimes varm fails for unknown reason
            means_per_cluster = adata_sc.var
        means_per_cluster = means_per_cluster[
            [
                f"means_per_cluster_mu_fg_{i}"
                for i in adata_sc.uns["mod"]["factor_names"]
            ]
        ].copy()
        means_per_cluster.columns = adata_sc.uns["mod"]["factor_names"]
    else:
        means_per_cluster = compute_cluster_averages(
            adata_sc,
            labels="label",
            layer=None,
            use_raw=False,
        )

    # SPATIAL MAPPING
    # find shared genes and subset both anndata and reference signatures
    intersect = np.intersect1d(adata.var_names, means_per_cluster.index)
    adata = adata[:, intersect].copy()
    means_per_cluster = means_per_cluster.loc[intersect, :].copy()

    # prepare anndata for cell2location model
    adata.obs["sample"] = "all"
    Cell2location.setup_anndata(adata=adata, batch_key="sample")
    cell2location_kwargs = dict(
        cell_state_df=means_per_cluster,
        # the expected average cell abundance: tissue-dependent
        # hyper-prior which can be estimated from paired histology:
        # here = average in the simulated dataset
        N_cells_per_location=n_cells_per_location,
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection:
        detection_alpha=detection_alpha,
    )
    if amortised:
        cell2location_kwargs["amortised"] = True
        cell2location_kwargs["encoder_mode"] = "multiple"
        cell2location_kwargs["encoder_kwargs"] = {
            "dropout_rate": 0.1,
            "n_hidden": {
                "single": 256,
                "n_s_cells_per_location": 10,
                "b_s_groups_per_location": 10,
                "z_sr_groups_factors": 64,
                "w_sf": 256,
                "detection_y_s": 20,
            },
            "use_batch_norm": False,
            "use_layer_norm": True,
            "n_layers": 1,
            "activation_fn": ELU,
        }
    # create and train the model
    st_model = Cell2location(adata, **cell2location_kwargs)
    st_model.train(
        max_epochs=max_epochs_st,
        # train using full data (batch_size=None)
        batch_size=st_batch_size,
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
            "batch_size": st_batch_size,
        },
    )

    adata.obsm["proportions_pred"] = adata.obsm["q05_cell_abundance_w_sf"].values
    adata.obsm["proportions_pred"] /= adata.obsm["proportions_pred"].sum(axis=1)[
        :, None
    ]
    adata.uns["method_code_version"] = check_version("cell2location")
    return adata


@_cell2location_method(
    method_name="Cell2location (detection_alpha=20, reference hard-coded)"
)
def cell2location_detection_alpha_20(
    adata,
    detection_alpha=20,
    n_cells_per_location=20,
    hard_coded_reference=True,
    amortised=False,
    num_samples=None,
    sc_batch_size=2500,
    st_batch_size=None,
    test: bool = False,
    max_epochs_sc: Optional[int] = None,
    max_epochs_st: Optional[int] = None,
):
    return _cell2location(
        adata,
        detection_alpha=detection_alpha,
        n_cells_per_location=n_cells_per_location,
        hard_coded_reference=hard_coded_reference,
        amortised=amortised,
        num_samples=num_samples,
        sc_batch_size=sc_batch_size,
        st_batch_size=st_batch_size,
        test=test,
        max_epochs_sc=max_epochs_sc,
        max_epochs_st=max_epochs_st,
    )


@_cell2location_method(
    method_name="Cell2location (detection_alpha=1, reference hard-coded)"
)
def cell2location_detection_alpha_1(
    adata,
    detection_alpha=1,
    n_cells_per_location=20,
    hard_coded_reference=True,
    amortised=False,
    num_samples=None,
    sc_batch_size=2500,
    st_batch_size=None,
    test: bool = False,
    max_epochs_sc: Optional[int] = None,
    max_epochs_st: Optional[int] = None,
):
    return _cell2location(
        adata,
        detection_alpha=detection_alpha,
        n_cells_per_location=n_cells_per_location,
        hard_coded_reference=hard_coded_reference,
        amortised=amortised,
        num_samples=num_samples,
        sc_batch_size=sc_batch_size,
        st_batch_size=st_batch_size,
        test=test,
        max_epochs_sc=max_epochs_sc,
        max_epochs_st=max_epochs_st,
    )


@_cell2location_method(
    method_name="Cell2location (detection_alpha=20, reference NB without batch info)"
)
def cell2location_detection_alpha_20_nb(
    adata,
    detection_alpha=20,
    n_cells_per_location=20,
    hard_coded_reference=False,
    amortised=False,
    num_samples=None,
    sc_batch_size=2500,
    st_batch_size=None,
    test: bool = False,
    max_epochs_sc: Optional[int] = None,
    max_epochs_st: Optional[int] = None,
):
    return _cell2location(
        adata,
        detection_alpha=detection_alpha,
        n_cells_per_location=n_cells_per_location,
        hard_coded_reference=hard_coded_reference,
        amortised=amortised,
        num_samples=num_samples,
        sc_batch_size=sc_batch_size,
        st_batch_size=st_batch_size,
        test=test,
        max_epochs_sc=max_epochs_sc,
        max_epochs_st=max_epochs_st,
    )


@_cell2location_method(
    method_name="Cell2location (detection_alpha=200, reference hard-coded)"
)
def cell2location_detection_alpha_200(
    adata,
    detection_alpha=200,
    n_cells_per_location=20,
    hard_coded_reference=True,
    amortised=False,
    num_samples=None,
    sc_batch_size=2500,
    st_batch_size=None,
    test: bool = False,
    max_epochs_sc: Optional[int] = None,
    max_epochs_st: Optional[int] = None,
):
    return _cell2location(
        adata,
        detection_alpha=detection_alpha,
        n_cells_per_location=n_cells_per_location,
        hard_coded_reference=hard_coded_reference,
        amortised=amortised,
        num_samples=num_samples,
        sc_batch_size=sc_batch_size,
        st_batch_size=st_batch_size,
        test=test,
        max_epochs_sc=max_epochs_sc,
        max_epochs_st=max_epochs_st,
    )


@_cell2location_method(
    method_name="Cell2location, amortised (detection_alpha=20, reference hard-coded)"
)
def cell2location_amortised_detection_alpha_20(
    adata,
    detection_alpha=20,
    n_cells_per_location=20,
    hard_coded_reference=True,
    amortised=True,
    num_samples=None,
    sc_batch_size=2500,
    st_batch_size=1024,
    test: bool = False,
    max_epochs_sc: Optional[int] = None,
    max_epochs_st: Optional[int] = None,
):
    return _cell2location(
        adata,
        detection_alpha=detection_alpha,
        n_cells_per_location=n_cells_per_location,
        hard_coded_reference=hard_coded_reference,
        amortised=amortised,
        num_samples=num_samples,
        sc_batch_size=sc_batch_size,
        st_batch_size=st_batch_size,
        test=test,
        max_epochs_sc=max_epochs_sc,
        max_epochs_st=max_epochs_st,
    )
