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

def download_file(url, destination, max_retries=5):
    """Download a file from a URL to a destination with progress bar and resume capability."""
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(destination), exist_ok=True)
    
    # Get file size from server
    response = requests.head(url)
    total_size_server = int(response.headers.get('content-length', 0))
    
    # Check if file exists and is complete
    if os.path.exists(destination):
        file_size = os.path.getsize(destination)
        if file_size == total_size_server:
            logger.info(f"File already exists and is complete at {destination}, skipping download.")
            return
        elif file_size < total_size_server:
            logger.info(f"Resuming download from byte {file_size} of {total_size_server}")
        else:
            logger.warning(f"Existing file is larger than expected. Restarting download.")
            file_size = 0
    else:
        file_size = 0
    
    # Try to download with retries
    for attempt in range(max_retries):
        try:
            # Set up headers for resume
            headers = {}
            mode = 'wb'
            if file_size > 0:
                headers['Range'] = f'bytes={file_size}-'
                mode = 'ab'  # Append to existing file
            
            # Make request with resume header if applicable
            response = requests.get(url, headers=headers, stream=True, timeout=60)
            
            # Handle resume response or normal response
            if response.status_code == 206:  # Partial content
                content_length = int(response.headers.get('content-length', 0))
                total_size = content_length + file_size
            elif response.status_code == 200:  # OK
                total_size = int(response.headers.get('content-length', 0))
                file_size = 0  # Reset file size for progress tracking
            else:
                logger.error(f"Unexpected status code: {response.status_code}")
                continue
            
            # Download with progress bar
            with open(destination, mode) as f:
                with tqdm(
                    total=total_size,
                    initial=file_size,
                    unit='iB',
                    unit_scale=True,
                    unit_divisor=1024,
                    desc=destination
                ) as bar:
                    for chunk in response.iter_content(chunk_size=1024*1024):  # 1MB chunks
                        if chunk:
                            size = f.write(chunk)
                            bar.update(size)
            
            # Verify file size
            if os.path.getsize(destination) >= total_size:
                logger.info(f"Download completed: {destination}")
                return
            else:
                logger.warning(f"Downloaded file size mismatch. Retrying...")
                file_size = os.path.getsize(destination)
                
        except (requests.exceptions.RequestException, IOError) as e:
            logger.warning(f"Download error (attempt {attempt+1}/{max_retries}): {str(e)}")
            # Update file size for next attempt
            if os.path.exists(destination):
                file_size = os.path.getsize(destination)
        
        # Wait before retrying
        if attempt < max_retries - 1:
            wait_time = 2 ** attempt  # Exponential backoff
            logger.info(f"Waiting {wait_time} seconds before retrying...")
            import time
            time.sleep(wait_time)
    
    raise Exception(f"Failed to download {url} after {max_retries} attempts")

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

def filter_op3_data(adata):
    """
    Filter the OP3 dataset based on specific criteria for each small molecule and cell type.
    """
    logger.info("Applying OP3-specific filtering criteria")
    
    # Create a boolean mask for filtering observations
    obs_filt = np.ones(adata.n_obs, dtype=bool)
    
    # Alvocidib only T cells in only 2 donors, remove
    obs_filt = obs_filt & (adata.obs['sm_name'] != "Alvocidib")
    
    # BMS-387032 - one donor with only T cells, two other consistent, but only 2 cell types
    # Leave the 2 cell types in, remove donor 2 with only T cells
    obs_filt = obs_filt & ~((adata.obs['sm_name'] == "BMS-387032") & (adata.obs['donor_id'] == "Donor 2"))
    
    # BMS-387032 remove myeloid cells and B cells
    obs_filt = obs_filt & ~((adata.obs['sm_name'] == "BMS-387032") & 
                           adata.obs['cell_type'].isin(["B cells", "Myeloid cells"]))
    
    # CGP 60474 has only T cells left, remove
    obs_filt = obs_filt & (adata.obs['sm_name'] != "CGP 60474")
    
    # Canertinib - the variation of Myeloid cell proportions is very large, skip Myeloid
    obs_filt = obs_filt & ~((adata.obs['sm_name'] == "Canertinib") & 
                           (adata.obs['cell_type'] == "Myeloid cells"))
    
    # Foretinib - large variation in Myeloid cell proportions (some in T cells), skip Myeloid
    obs_filt = obs_filt & ~((adata.obs['sm_name'] == "Foretinib") & 
                           (adata.obs['cell_type'] == "Myeloid cells"))
    
    # Ganetespib (STA-9090) - donor 2 has no Myeloid and small NK cells proportions
    # Skip Myeloid, remove donor 2
    obs_filt = obs_filt & ~((adata.obs['sm_name'] == "Ganetespib (STA-9090)") & 
                           (adata.obs['donor_id'] == "Donor 2"))
    
    # IN1451 - donor 2 has no NK or B, remove Donor 2
    obs_filt = obs_filt & ~((adata.obs['sm_name'] == "IN1451") & 
                           (adata.obs['donor_id'] == "Donor 2"))
    
    # Navitoclax - donor 3 doesn't have B cells and has different T and Myeloid proportions
    # Remove donor 3
    obs_filt = obs_filt & ~((adata.obs['sm_name'] == "Navitoclax") & 
                           (adata.obs['donor_id'] == "Donor 3"))
    
    # PF-04691502 remove Myeloid (only present in donor 3)
    obs_filt = obs_filt & ~((adata.obs['sm_name'] == "PF-04691502") & 
                           (adata.obs['cell_type'] == "Myeloid cells"))
    
    # Proscillaridin A;Proscillaridin-A remove Myeloid, since the variation is very high (4x)
    obs_filt = obs_filt & ~((adata.obs['sm_name'] == "Proscillaridin A;Proscillaridin-A") & 
                           (adata.obs['cell_type'] == "Myeloid cells"))
    
    # R428 - skip NK due to high variation (close to 3x)
    obs_filt = obs_filt & ~((adata.obs['sm_name'] == "R428") & 
                           (adata.obs['cell_type'] == "NK cells"))
    
    # UNII-BXU45ZH6LI - remove due to large variation across all cell types and missing cell types
    obs_filt = obs_filt & (adata.obs['sm_name'] != "UNII-BXU45ZH6LI")
    
    # Apply the filter
    filtered_adata = adata[obs_filt].copy()
    
    logger.info(f"Filtered data from {adata.n_obs} to {filtered_adata.n_obs} cells")
    
    return filtered_adata

def filter_by_counts(adata, par):
    """Filter cells and genes by count thresholds."""
    logger.info("Filtering cells and genes by count thresholds")
    n_cells_before, n_genes_before = adata.shape
    
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
    
    # Apply OP3-specific filtering
    adata = filter_op3_data(adata)
    
    # Filter by parameters
    if par["donor_id"] is not None:
        logger.info(f"Filtering for donor_id: {par['donor_id']}")
        adata = adata[adata.obs["donor_id"] == par["donor_id"]]
    
    if par["cell_type"] is not None:
        logger.info(f"Filtering for cell_type: {par['cell_type']}")
        adata = adata[adata.obs["cell_type"] == par["cell_type"]]
    
    if par["perturbation"] is not None:
        logger.info(f"Filtering for perturbation: {par['perturbation']}")
        adata = adata[adata.obs["perturbation"] == par["perturbation"]]
    
    # Filter cells and genes
    adata = filter_by_counts(adata, par)
    
    # Move X to layers and normalize
    move_x_to_layers(adata)
    
    # Add dataset metadata
    add_metadata_to_uns(adata, par)
    
    # Print summary and save
    print_summary(adata)
    write_anndata(adata, par)
    
    logger.info("Done")

if __name__ == "__main__":
    main(par, meta)