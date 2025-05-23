import os
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import requests
import logging
from urllib.parse import urlparse

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

## VIASH START
par = {
    "input": "resources/neurips-2023-raw/sc_counts.h5ad",
    "var_feature_name": "index",
    "data_type": "sc",
    "donor_id": None,
    "cell_type": None,
    "perturbation": None,
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
    """Download a file from a URL to a destination."""
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
            
            logger.info(f"Downloading {url} to {destination}")
            
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
            
            with open(destination, mode) as f:
                for chunk in response.iter_content(chunk_size=1024*1024):  # 1MB chunks
                    if chunk:
                        f.write(chunk)
            
            # Normal verification code
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

def move_x_to_layers(adata):
    """Move .X to .layers['counts'] and set X to None."""
    logger.info("Moving .X to .layers['counts']")
    adata.layers["counts"] = adata.X.copy()
    adata.X = None

def add_metadata_to_uns(adata, par):
    """Add standardized metadata to .uns."""
    logger.info("Adding metadata to .uns")
    adata.uns['dataset_id'] = par["dataset_id"]
    adata.uns['dataset_name'] = par["dataset_name"]
    adata.uns['dataset_summary'] = par["dataset_summary"]
    adata.uns['dataset_description'] = par["dataset_description"]
    adata.uns['dataset_organism'] = 'Homo sapiens'
    adata.uns['dataset_url'] = par["dataset_url"]
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

# Instead of defining main() and calling it at the end, write the code directly
logger.info("Starting OP3 loader")

# Load the data
logger.info(f"Loading data from {par['input']}")

if os.path.isfile(par["input"]):
    try:
        adata = sc.read_h5ad(par["input"])
    except:
        raise ValueError('The downloaded dataset does not match the .h5ad type. ')
elif (parsed := urlparse(par["dataset_url"])) and (parsed.scheme and parsed.netloc):
    download_file(par["dataset_url"], par["input"], max_retries=5)
    try:
        adata = sc.read_h5ad(par["input"])
    except:
        raise ValueError('The downloaded dataset does not match the .h5ad type. ')
else:
    raise ValueError('The dataset is not found. Please enter another link/path to it.')

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


# Move X to layers and normalize
move_x_to_layers(adata)

# Add dataset metadata
add_metadata_to_uns(adata, par)

# Add a feature name
logger.info("Setting .var['feature_name']")

if par["var_feature_name"] == "index":
    adata.var["feature_name"] = adata.var.index
else:
    if par["var_feature_name"] in adata.var:
        adata.var["feature_name"] = adata.var[par["feature_name"]]
        del adata.var[par["feature_name"]]
    else:
        logger.info(f"Warning: key '{par['var_feature_name']}' could not be found in adata.var.")

# Print summary and save
print_summary(adata)
write_anndata(adata, par)

logger.info("Done")
