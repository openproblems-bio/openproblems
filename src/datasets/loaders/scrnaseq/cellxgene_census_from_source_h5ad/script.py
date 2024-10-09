import sys
import cellxgene_census
import scanpy as sc
import tempfile

## VIASH START
par = {
    "input_id": "0895c838-e550-48a3-a777-dbcd35d30272",
    "obs_batch": [ "donor_id" ],
    "obs_batch_separator": "+",
    "dataset_name": "pretty name",
    "dataset_url": "url",
    "dataset_reference": "ref",
    "dataset_summary": "summ",
    "dataset_description": "desc",
    "dataset_organism": "mus_musculus",
    "output": "output.h5ad",
    "output_compression": "gzip",
}
meta = {"resources_dir": "src/common/helper_functions"}
## VIASH END

sys.path.append(meta["resources_dir"])

from setup_logger import setup_logger
logger = setup_logger()

def get_anndata(par):
    with tempfile.TemporaryDirectory() as tmp:
        path = tmp + "/source.h5ad"
        logger.info("Downloading source h5ad for dataset '%s' to '%s'.", par["input_id"], path)
        cellxgene_census.download_source_h5ad(par["input_id"], path)
        return sc.read_h5ad(path)

def filter_by_counts(adata, par):
    logger.info("Remove cells with few counts and genes with few counts.")
    t0 = adata.shape
    # remove cells with few counts and genes with few counts
    if par["cell_filter_min_counts"]:
        sc.pp.filter_cells(adata, min_counts=par["cell_filter_min_counts"])
    if par["cell_filter_min_genes"]:
        sc.pp.filter_cells(adata, min_genes=par["cell_filter_min_genes"])
    if par["gene_filter_min_counts"]:
        sc.pp.filter_genes(adata, min_counts=par["gene_filter_min_counts"])
    if par["gene_filter_min_cells"]:
        sc.pp.filter_genes(adata, min_cells=par["gene_filter_min_cells"])
    t1 = adata.shape
    logger.info("Removed %s cells and %s genes.", (t0[0] - t1[0]), (t0[1] - t1[1]))

def move_x_to_layers(adata):
    logger.info("Move .X to .layers['counts']")
    adata.layers["counts"] = adata.X
    adata.X = None

def add_batch_to_obs(adata, par):
    logger.info("Add batch to the AnnData object.")
    if par["obs_batch"]:
        # fetch batch columns from obs
        cols = [adata.obs[key] for key in par["obs_batch"]]
        
        # join cols
        obs_batch = [par["obs_batch_separator"].join(row) for row in zip(*cols)]

        # store in adata
        adata.obs["batch"] = obs_batch

def add_metadata_to_uns(adata, par):
    logger.info("Add metadata to the AnnData object.")
    for key in ["dataset_id", "dataset_name", "dataset_url", "dataset_reference", "dataset_summary", "dataset_description", "dataset_organism"]:
        adata.uns[key] = par[key]

def print_unique(adata, column):
    if column not in adata.obs.columns:
        logger.info(f"Column {column} not found in obs")
        return
    formatted = "', '".join(adata.obs[column].unique())
    logger.info(f"Unique {column}: ['{formatted}']")

def print_summary(adata):
    logger.info(f"Resulting dataset: {adata}")

    logger.info("Summary of dataset:")
    print_unique(adata, "assay")
    print_unique(adata, "assay_ontology_term_id")
    print_unique(adata, "cell_type")
    print_unique(adata, "cell_type_ontology_term_id")
    print_unique(adata, "dataset_id")
    print_unique(adata, "development_stage")
    print_unique(adata, "development_stage_ontology_term_id")
    print_unique(adata, "disease")
    print_unique(adata, "disease_ontology_term_id")
    print_unique(adata, "tissue")
    print_unique(adata, "tissue_ontology_term_id")
    print_unique(adata, "tissue_general")
    print_unique(adata, "tissue_general_ontology_term_id")

def write_anndata(adata, par):
    logger.info("Writing AnnData object to '%s'", par["output"])

    adata.write_h5ad(par["output"], compression=par["output_compression"])

def main(par, meta):
    adata = get_anndata(par)

    logger.info("AnnData: %s", str(adata))

    # remove cells with few counts and genes with few counts
    filter_by_counts(adata, par)

    # this is not needed in source h5ads
    # # use feature_id as var_names
    # adata.var_names = adata.var["feature_id"]

    # move .X to .layers["counts"]
    move_x_to_layers(adata)

    # add batch to obs
    add_batch_to_obs(adata, par)

    # add metadata to uns
    add_metadata_to_uns(adata, par)

    # print summary
    print_summary(adata)

    # write output to file
    write_anndata(adata, par)


if __name__ == "__main__":
    main(par, meta)
