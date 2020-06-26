import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import scprep
import tempfile
import os


def load_scicar(
    rna_url,
    rna_cells_url,
    rna_genes_url,
    atac_url,
    atac_cells_url,
    atac_genes_url,
    test=False,
):
    rna_cells = pd.read_csv(rna_cells_url, low_memory=False)["sample"]
    rna_genes = pd.read_csv(rna_genes_url, low_memory=False)["gene_id"]
    atac_cells = pd.read_csv(atac_cells_url, low_memory=False)["sample"]
    atac_genes = pd.read_csv(atac_genes_url, low_memory=False, index_col=0)["peak"]

    with tempfile.TemporaryDirectory() as tempdir:
        rna_file = os.path.join(tempdir, "rna.mtx.gz")
        scprep.io.download.download_url(rna_url, rna_file)
        rna_data = scprep.io.load_mtx(rna_file, cell_axis="col").tocsr()
        atac_file = os.path.join(tempdir, "atac.mtx.gz")
        scprep.io.download.download_url(atac_url, atac_file)
        atac_data = scprep.io.load_mtx(atac_file, cell_axis="col").tocsr()

    if test:
        remove_genes = rna_data.sum(axis=0).A.flatten() < 1
        rna_genes, rna_data = rna_genes[~remove_genes], rna_data[:, ~remove_genes]
        remove_genes = atac_data.sum(axis=0).A.flatten() < 1
        atac_genes, atac_data = atac_genes[~remove_genes], atac_data[:, ~remove_genes]

        rna_genes, rna_data = rna_genes[:200], rna_data[:, :200]
        atac_genes, atac_data = atac_genes[:400], atac_data[:, :400]

    rna_data, rna_cells = scprep.filter.filter_empty_cells(rna_data, rna_cells)
    atac_data, atac_cells = scprep.filter.filter_empty_cells(atac_data, atac_cells)

    common_cells = np.intersect1d(rna_cells, atac_cells)
    if test:
        common_cells = common_cells[:100]

    rna_subset = np.isin(rna_cells, common_cells)
    rna_cells, rna_data = rna_cells[rna_subset], rna_data[rna_subset]
    rna_order = np.argsort(rna_cells.to_numpy())
    rna_cells, rna_data = rna_cells.iloc[rna_order], rna_data[rna_order]

    atac_subset = np.isin(atac_cells, common_cells)
    atac_cells, atac_data = atac_cells[atac_subset], atac_data[atac_subset]
    atac_order = np.argsort(atac_cells.to_numpy())
    atac_cells, atac_data = atac_cells.iloc[atac_order], atac_data[atac_order]

    adata = anndata.AnnData(
        rna_data, obs=pd.DataFrame(index=rna_cells), var=pd.DataFrame(index=rna_genes),
    )
    adata.uns["mode2"] = anndata.AnnData(
        atac_data,
        obs=pd.DataFrame(index=atac_cells),
        var=pd.DataFrame(index=atac_genes),
    )

    adata.raw = adata
    adata.uns["mode2"].raw = adata.uns["mode2"]

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.normalize_total(adata.uns["mode2"], target_sum=1e4)

    sc.pp.log1p(adata)
    sc.pp.log1p(adata.uns["mode2"])
    return adata
