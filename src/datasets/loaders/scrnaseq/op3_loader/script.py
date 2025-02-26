import os
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import requests
from tqdm import tqdm # liberary for displaying progress bars
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

## VIASH START
par = {
    "data_type": "sc",
    "donor_id": None,
    "cell_type": None,
    "perturbation": None,
    "min_cells": 3,
    "min_genes": 200,
    "dataset_id": "op3",
    "dataset_name": "OP3: single-cell multimodal dataset in PBMCs for perturbation prediction benchmarking",
    "dataset_summary": "The Open Problems Perurbation Prediction (OP3) dataset with small molecule perturbations in PBMCs",
    "dataset_description": "The OP3 dataset is to-date the largest single-cell small molecule perturbation dataset in primary tissue with multiple donor replicates.",
    "output": "output.h5ad",
    "output_compression": "gzip"
}
meta = {}
## VIASH END

def download_file(url, destination):
    """Download a file from a URL to a destination with progress bar."""
    if os.path.exists(destination):
        logger.info(f"File already exists at {destination}, skipping download.")
        return
    
    logger.info(f"Downloading {url} to {destination}")
    response = requests.get(url, stream=True)
    response.raise_for_status()
    
    total_size = int(response.headers.get('content-length', 0))
    block_size = 1024  # 1 Kibibyte
    
    with open(destination, 'wb') as file, tqdm(
            desc=destination,
            total=total_size,
            unit='iB',
            unit_scale=True,
            unit_divisor=1024,
        ) as bar:
        for data in response.iter_content(block_size):
            size = file.write(data)
            bar.update(size)

def download_op3_data():
    """Download the OP3 dataset from GEO."""
    base_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279945/suppl/"
    filename = "GSE279945_sc_counts_processed.h5ad"
    
    url = f"{base_url}{filename}"
    cache_dir = os.path.join(os.path.expanduser("~"), ".cache", "op3_loader")
    os.makedirs(cache_dir, exist_ok=True)
    
    destination = os.path.join(cache_dir, filename)
    download_file(url, destination)
    
    return destination

def filter_adata(adata, donor_id=None, cell_type=None, perturbation=None):
    """Filter AnnData object based on donor_id, cell_type, and perturbation."""
    n_cells_before = adata.n_obs
    
    if donor_id is not None:
        logger.info(f"Filtering for donor_id: {donor_id}")
        adata = adata[adata.obs['donor_id'] == donor_id]
        
    if cell_type is not None:
        logger.info(f"Filtering for cell_type: {cell_type}")
        adata = adata[adata.obs['cell_type'] == cell_type]
        
    if perturbation is not None:
        logger.info(f"Filtering for perturbation: {perturbation}")
        adata = adata[adata.obs['perturbation'] == perturbation]
    
    n_cells_after = adata.n_obs
    logger.info(f"Filtered from {n_cells_before} to {n_cells_after} cells")
        
    return adata

def filter_by_counts(adata, par):
    """Filter cells and genes by count thresholds."""
    logger.info("Filtering cells and genes by count thresholds")
    n_cells_before, n_genes_before = adata.shape
    
    # Basic filtering
    sc.pp.filter_cells(adata, min_genes=par["min_genes"])
    sc.pp.filter_genes(adata, min_cells=par["min_cells"])
    
    n_cells_after, n_genes_after = adata.shape
    logger.info(f"Removed {n_cells_before - n_cells_after} cells and {n_genes_before - n_genes_after} genes")
    
    return adata

def move_x_to_layers(adata):
    """Move .X to .layers['counts'] and normalize data."""
    logger.info("Moving .X to .layers['counts']")
    adata.layers["counts"] = adata.X.copy()
    
    # Normalize and log transform
    logger.info("Normalizing and log-transforming data")
    # Normalization scales gene expression values to be comparable between cells
    # by accounting for differences in sequencing depth
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

def add_metadata_to_uns(adata, par):
    """Add standardized metadata to .uns."""
    logger.info("Adding metadata to .uns")
    adata.uns['dataset_id'] = par["dataset_id"]
    adata.uns['dataset_name'] = par["dataset_name"]
    adata.uns['dataset_summary'] = par["dataset_summary"]
    adata.uns['dataset_description'] = par["dataset_description"]
    adata.uns['dataset_organism'] = 'Homo sapiens'
    adata.uns['dataset_url'] = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279945'
    adata.uns['dataset_reference'] = 'GSE279945'
    adata.uns['dataset_version'] = '1.0.0'
    adata.uns['processing_status'] = 'processed'

def print_unique(adata, column):
    """Print unique values in a column."""
    if column in adata.obs:
        values = adata.obs[column].unique()
        if len(values) <= 10:
            formatted = "', '".join(values)
            logger.info(f"Unique {column}: ['{formatted}']")
        else:
            logger.info(f"Unique {column}: {len(values)} values")

def print_summary(adata):
    """Print a summary of the dataset."""
    logger.info(f"Dataset shape: {adata.shape}")
    logger.info(f"Number of cells: {adata.n_obs}")
    logger.info(f"Number of genes: {adata.n_vars}")
    
    # Print unique values for key columns
    for column in ['donor_id', 'cell_type', 'perturbation']:
        print_unique(adata, column)
        
    logger.info(f"Layers: {list(adata.layers.keys())}")
    logger.info(f"Metadata: {list(adata.uns.keys())}")

def write_anndata(adata, par):
    """Write AnnData object to file."""
    logger.info(f"Writing AnnData object to '{par['output']}'")
    adata.write_h5ad(par["output"], compression=par["output_compression"])

def main(par, meta):
    """Main function."""
    logger.info("Starting OP3 loader")
    
    # Download the data
    data_path = download_op3_data()
    
    # Load the data
    logger.info(f"Loading data from {data_path}")
    adata = sc.read_h5ad(data_path)
    
    # Filter by donor_id if specified
    if par["donor_id"] is not None:
        logger.info(f"Filtering for donor_id: {par['donor_id']}")
        adata = adata[adata.obs["donor_id"] == par["donor_id"]]
    
    # Filter by cell_type if specified
    if par["cell_type"] is not None:
        logger.info(f"Filtering for cell_type: {par['cell_type']}")
        adata = adata[adata.obs["cell_type"] == par["cell_type"]]
    
    # Filter by perturbation if specified
    if par["perturbation"] is not None:
        logger.info(f"Filtering for perturbation: {par['perturbation']}")
        adata = adata[adata.obs["perturbation"] == par["perturbation"]]
    
    # Filter cells and genes
    logger.info("Filtering cells and genes")
    sc.pp.filter_cells(adata, min_genes=par["min_genes"])
    sc.pp.filter_genes(adata, min_cells=par["min_cells"])
    
    # Move X to layers and normalize
    move_x_to_layers(adata)
    
    # Add dataset metadata
    logger.info("Adding dataset metadata")
    adata.uns["dataset_id"] = par["dataset_id"]
    adata.uns["dataset_name"] = par["dataset_name"]
    adata.uns["dataset_summary"] = par["dataset_summary"]
    adata.uns["dataset_description"] = par["dataset_description"]
    
    # Write output
    logger.info(f"Writing output to {par['output']}")
    adata.write_h5ad(par["output"], compression=par["output_compression"])
    
    logger.info("Done")

if __name__ == "__main__":
    main(par, meta)