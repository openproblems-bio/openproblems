from . import utils

import os
import requests
import scanpy as sc
import scprep
import tempfile

COLLECTION_ID = "0b9d8a04-bb9d-44da-aa27-705bb65b54eb"
DOMAIN = "cellxgene.cziscience.com"
API_BASE = f"https://api.{DOMAIN}"
METHOD_ALIASES = {
    "10x 3' v2": "droplet",
    "Smart-seq2": "facs"
}


def matching_dataset(dataset, method_list, organ_list):
    # we want to skip combined datasets
    if len(dataset["tissue"]) > 1:
        return False
    if dataset["tissue"][0]["label"] not in organ_list:
        return False
    method = METHOD_ALIASES[dataset["assay"][0]["label"]]
    return method in method_list


def load_raw_counts(dataset):
    dataset_id = dataset["id"]
    assets_path = f"/curation/v1/collections/{COLLECTION_ID}/datasets/{dataset_id}/assets"
    url = f"{API_BASE}{assets_path}"
    res = requests.get(url=url)
    assets = res.json()
    assets = [asset for asset in assets if asset["filetype"] == "H5AD"]
    asset = assets[0]

    filename = f"{COLLECTION_ID}_{dataset_id}_{asset['filename']}"
    with tempfile.TemporaryDirectory() as tempdir:
        filepath = os.path.join(tempdir, filename)
        scprep.io.download.download_url(asset["presigned_url"], filepath)
        adata = sc.read_h5ad(filepath)

    utils.filter_genes_cells(adata)
    # raw counts are in `raw`
    return adata.raw.to_adata()


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

    if not method_list:
        raise ValueError("`method_list' argument not provided")
    unknown_methods = set(method_list) - set(["facs", "droplet"])
    if unknown_methods:
        raise ValueError(
            f"Unknown methods provided in `method_list': {','.join(unknown_methods)}. "
            "Known methods are `facs' and `droplet'"
        )
    organ_list = [x.lower() for x in organ_list]

    datasets_path = f"/curation/v1/collections/{COLLECTION_ID}"
    url = f"{API_BASE}{datasets_path}"
    res = requests.get(url=url)
    datasets = res.json()["datasets"]

    adata_list = []
    for dataset in datasets:
        if matching_dataset(dataset, method_list, organ_list):
            adata_list.append(load_raw_counts(dataset))

    if len(adata_list) > 1:
        adata = adata_list[0].concatenate(adata_list[1:], join="outer")
    else:
        adata = adata_list[0]

    if test:
        sc.pp.subsample(adata, n_obs=500)
        adata = adata[:, :1000]
        utils.filter_genes_cells(adata)
    return adata    
