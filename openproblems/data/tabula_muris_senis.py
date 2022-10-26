from . import utils

import anndata as ad
import os
import requests
import scanpy as sc
import scprep
import tempfile

COLLECTION_ID = "0b9d8a04-bb9d-44da-aa27-705bb65b54eb"
DOMAIN = "cellxgene.cziscience.com"
API_BASE = f"https://api.{DOMAIN}"
METHOD_ALIASES = {"10x 3' v2": "droplet", "Smart-seq2": "facs"}


def check_unknown_organs(datasets, organ_list):
    known_organs = set([t["label"] for d in datasets for t in d["tissue"]])
    unknown_organs = set(organ_list) - known_organs
    if unknown_organs:
        raise ValueError(
            f"Unknown organs provided in `organ_list': {', '.join(unknown_organs)}. "
            f"Known organs are {', '.join(known_organs)}"
        )


def matching_dataset(dataset, method_list, organ_list):
    # if dataset has multiple methods, skip it
    if len(dataset["assay"]) > 1:
        return False

    # if dataset has multiple tissues, skip it
    if len(dataset["tissue"]) > 1:
        return False

    method = dataset["assay"][0]["label"]
    method = METHOD_ALIASES[method]

    # if organ_list is not empty, check for specific tissue
    if len(organ_list) > 0 and dataset["tissue"][0]["label"] not in organ_list:
        return False

    # if method_list is not empty, check for specific method
    if len(method_list) > 0 and method not in method_list:
        return False

    return True


def load_raw_counts(dataset):
    dataset_id = dataset["id"]
    assets_path = (
        f"/curation/v1/collections/{COLLECTION_ID}/datasets/{dataset_id}/assets"
    )
    url = f"{API_BASE}{assets_path}"
    res = requests.get(url=url)
    assets = res.json()
    assets = [asset for asset in assets if asset["filetype"] == "H5AD"]
    assert len(assets) == 1
    asset = assets[0]

    filename = f"{COLLECTION_ID}_{dataset_id}_{asset['filename']}"
    with tempfile.TemporaryDirectory() as tempdir:
        filepath = os.path.join(tempdir, filename)
        scprep.io.download.download_url(asset["presigned_url"], filepath)
        adata = sc.read_h5ad(filepath)

    utils.filter_genes_cells(adata)
    # If `raw` exists, raw counts are there
    if getattr(adata, "raw", None) is not None:
        return adata.raw.to_adata()
    return adata


@utils.loader(
    data_url="https://tabula-muris-senis.ds.czbiohub.org/",
    data_reference="https://doi.org/10.1038/s41586-020-2496-1",
)
def load_tabula_muris_senis(test=False, method_list=None, organ_list=None):
    """Load tubula_muris_senis datasets into 1 anndata object based on user input.

    Input which methods and organs to create anndata object from.
    Returns a single anndata object with specified methods and organs.
    EX: load_tabula_muris_senis(method_list = ['facs', 'droplet'],
    organ_list = ['Skin', 'Fat']), and returns anndata for facs-skin, droplet-skin,
    and droplet-fat anndata sets. (no facs-fat dataset available)
    """

    if method_list is None:
        method_list = []
    if organ_list is None:
        organ_list = []
    method_list = [x.lower() for x in method_list]
    organ_list = [x.lower() for x in organ_list]

    unknown_methods = set(method_list) - set(["facs", "droplet"])
    if unknown_methods:
        raise ValueError(
            f"Unknown methods provided in `method_list': {','.join(unknown_methods)}. "
            "Known methods are `facs' and `droplet'"
        )

    datasets_path = f"/curation/v1/collections/{COLLECTION_ID}"
    url = f"{API_BASE}{datasets_path}"
    res = requests.get(url=url)
    datasets = res.json()["datasets"]
    check_unknown_organs(datasets, organ_list)

    adata_list = []
    for dataset in datasets:
        if matching_dataset(dataset, method_list, organ_list):
            adata_list.append(load_raw_counts(dataset))

    assert len(adata_list) > 0
    adata = ad.concat(adata_list, join="outer")

    if test:
        adata = utils.subsample_even(adata, n_obs=500, even_obs="method")
    return adata
