import anndata as ad
import numpy as np
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

## VIASH START
par = {
    "input": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279945/suppl/GSE279945_sc_counts_processed.h5ad",
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

def filter_op3_data(adata):
    """
    Filter the OP3 dataset based on specific criteria for each small molecule and cell type.
    """
    logger.info("Applying OP3-specific filtering criteria")
    
    # Original filter code follows...
    
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

def fix_splits(adata):
    """Fix splits after reannotation."""
    logger.info("Fix splits after reannotation")
    adata.obs["cell_type_orig_updated"] = adata.obs["cell_type_orig"].apply(lambda x: "T cells" if x.startswith("T ") else x)
    adata.obs["sm_cell_type_orig"] = adata.obs["sm_name"].astype(str) + "_" + adata.obs["cell_type_orig_updated"].astype(str)
    mapping_to_split = adata.obs.groupby("sm_cell_type_orig")["split"].apply(lambda x: x.unique()[0]).to_dict()
    adata.obs["sm_cell_type"] = adata.obs["sm_name"].astype(str) + "_" + adata.obs["cell_type"].astype(str)
    adata.obs["split"] = adata.obs["sm_cell_type"].map(mapping_to_split)
    adata.obs['control'] = adata.obs['split'].eq("control")
    return adata
    

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
logger.info(f"Loading data at {par['input']}")
adata = ad.read_h5ad(par["input"])

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

# Fix splits after reannotation
fix_splits(adata)

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
