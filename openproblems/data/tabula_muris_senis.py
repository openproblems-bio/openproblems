from . import utils

import os
import pandas as pd
import scanpy as sc
import scprep
import tempfile


def get_filenames_and_urls(url_df, method_list=None, organ_list=None):
    """Takes in dataframe and returns corresponding filename(s) and url(s).

    Takes in dataframe (with sample information stored), a list of methods,
    and a list of organs.
    Returns filenames and figshare URLs associated with inputs.
    If method_list or organ_list are None, do not filter based on that argument.
    """
    subset_df = url_df.copy()
    # If method_list specified, filter based on methods in list.
    if method_list:
        subset_df = subset_df.loc[subset_df.method.isin(method_list)]
    # If organ_list specified, filter based on organs in list.
    if organ_list:
        subset_df = subset_df.loc[subset_df.organ.isin(organ_list)]

    return subset_df


def make_anndata_from_filename_and_url(filename, url, test=False):
    """Takes in filename and url pair. Returns corresponding anndata object."""
    with tempfile.TemporaryDirectory() as tempdir:
        filepath = os.path.join(tempdir, filename)
        scprep.io.download.download_url(url, filepath)
        adata = sc.read_h5ad(filepath)
        utils.filter_genes_cells(adata)

    if test:
        sc.pp.subsample(adata, n_obs=100)
        adata = adata[:, :1000]
        utils.filter_genes_cells(adata)

    return adata


def make_anndata_list(subset_df, test):
    """Makes anndata from filename/url pair. Adds to list of anndatas.

    Input dataframe that contains filenames and urls to make anndatas from.
    Returns a list of anndata objects.
    """
    adata_list = []
    for i in range(len(subset_df)):
        row = subset_df.iloc[i]
        adata_list.append(
            make_anndata_from_filename_and_url(row.filename, row.figshare_url)
        )
    if test:
        return adata_list[0]
    return adata_list


def combine_anndata(anndata_list):
    """Combines list of anndata objects and returns 1 anndata object.

    Takes in a list of anndata objects.
    Returns 1 anndata with all anndata objects combined.
    """

    # if anndata_list only contains 1 anndata object, will return that object
    anndata_concat = anndata_list[0].concatenate(anndata_list[1:])
    return anndata_concat


@utils.loader
def load_tabula_muris_senis(test=False, method_list=None, organ_list=None):
    """Load tubula_muris_senis datasets into 1 anndata object based on user input.

    Input which methods and organs to create anndata object from.
    Returns a single anndata object with specified methods and organs.
    EX: load_tabula_muris_senis(method_list = ['facs', 'droplet'],
    organ_list = ['Skin', 'Fat']), and returns anndata for facs-skin, droplet-skin,
    and droplet-fat anndata sets. (no facs-fat dataset available)
    """

    # df containing figshare links, method of collection, and organ for each
    # tabula muris dataset
    url_df = pd.read_csv(
        os.path.join(
            os.path.dirname(__file__),
            "tabula_muris_senis_data_objects",
            "tabula_muris_senis_data_objects.csv",
        ),
        header=0,
    )

    subset_df = get_filenames_and_urls(url_df, method_list, organ_list)
    adata_list = make_anndata_list(subset_df, test)
    adata_final = combine_anndata(adata_list)
    return adata_final
