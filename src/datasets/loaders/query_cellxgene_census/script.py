import os
import sys
import cellxgene_census

## VIASH START
par = {
    "input_uri": None,
    "census_version": "stable",
    "species": "homo_sapiens",
    "obs_value_filter": "is_primary_data == True and cell_type_ontology_term_id in ['CL:0000136', 'CL:1000311', 'CL:0002616'] and suspension_type == 'cell'",
    "cell_filter_grouping": ["dataset_id", "tissue", "assay", "disease", "cell_type"],
    "cell_filter_minimum_count": 100,
    "output": "output.h5ad",
    "output_compression": "gzip",
}
meta = {"resources_dir": "src/common/helper_functions"}
## VIASH END

sys.path.append(meta["resources_dir"])

from setup_logger import setup_logger
logger = setup_logger()

def connect_census(uri, census_version):
    """
    Connect to CellxGene Census or user-provided TileDBSoma object
    """
    ver = census_version or "stable"
    logger.info("Connecting to CellxGene Census at %s", f"'{uri}'" if uri else f"version '{ver}'")
    return cellxgene_census.open_soma(uri=uri, census_version=ver)

def get_anndata(census_connection, obs_value_filter, species):
    logger.info("Getting gene expression data based on %s query.", obs_value_filter)
    return cellxgene_census.get_anndata(
        census=census_connection, obs_value_filter=obs_value_filter, organism=species
    )


def add_cellcensus_metadata_obs(census_connection, query_data):
    logger.info("Adding extented metadata to gene expression data.")
    census_datasets = (
        census_connection["census_info"]["datasets"].read().concat().to_pandas()
    )

    query_data.obs.dataset_id = query_data.obs.dataset_id.astype("category")

    dataset_info = (
        census_datasets[
            census_datasets.dataset_id.isin(query_data.obs.dataset_id.cat.categories)
        ][
            [
                "collection_id",
                "collection_name",
                "collection_doi",
                "dataset_id",
                "dataset_title",
            ]
        ]
        .reset_index(drop=True)
        .apply(lambda x: x.astype("category"))
    )

    return query_data.obs.merge(dataset_info, on="dataset_id", how="left")


def cellcensus_cell_filter(query_data, cell_filter_grouping, cell_filter_minimum_count):
    t0 = query_data.shape
    query_data = query_data[
        query_data.obs.groupby(cell_filter_grouping)["soma_joinid"].transform("count")
        >= cell_filter_minimum_count
    ]
    t1 = query_data.shape
    logger.info(
        "Removed %s cells based on %s cell_filter_minimum_count of %s cell_filter_grouping."
        % ((t0[0] - t1[0]), cell_filter_minimum_count, cell_filter_grouping)
    )
    return query_data


def write_anndata(query_data, path, compression):
    logger.info("Writing AnnData object to '%s'", path)

    query_data.write_h5ad(path, compression=compression)

def print_unique(adata, column):
    formatted = "', '".join(adata.obs[column].unique())
    logger.info(f"Unique {column}: ['{formatted}']")

def print_summary(query_data):
    logger.info(f"Resulting dataset: {query_data}")

    logger.info("Summary of dataset:")
    print_unique(query_data, "assay")
    print_unique(query_data, "assay_ontology_term_id")
    print_unique(query_data, "cell_type")
    print_unique(query_data, "cell_type_ontology_term_id")
    print_unique(query_data, "dataset_id")
    print_unique(query_data, "development_stage")
    print_unique(query_data, "development_stage_ontology_term_id")
    print_unique(query_data, "disease")
    print_unique(query_data, "disease_ontology_term_id")
    print_unique(query_data, "tissue")
    print_unique(query_data, "tissue_ontology_term_id")
    print_unique(query_data, "tissue_general")
    print_unique(query_data, "tissue_general_ontology_term_id")

def main():
    # check arguments
    if (par["cell_filter_grouping"] is None) != (par["cell_filter_minimum_count"] is None):
        raise NotImplementedError(
            "You need to specify either both or none of the following parameters: cell_filter_grouping, cell_filter_minimum_count"
        )
    
    with connect_census(uri=par["input_uri"], census_version=par["census_version"]) as conn:
        query_data = get_anndata(conn, par["obs_value_filter"], par["species"])

        if par["add_collection_metadata"]:
            query_data.obs = add_cellcensus_metadata_obs(conn, query_data)

    if par["cell_filter_grouping"] is not None:
        query_data = cellcensus_cell_filter(
            query_data,
            par["cell_filter_grouping"],
            par["cell_filter_minimum_count"]
        )

    # use feature_id as var_names
    query_data.var_names = query_data.var["feature_id"]

    # print summary
    print_summary(query_data)

    # write output to file
    write_anndata(query_data, par["output"], par["output_compression"])


if __name__ == "__main__":
    main()
