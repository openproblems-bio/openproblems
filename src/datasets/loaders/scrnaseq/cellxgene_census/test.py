import sys
import os
import pytest
import anndata as ad
import numpy as np

## VIASH START
meta = {
    'resources_dir': './resources_test/',
    'executable': './target/docker/query/cellxgene_census',
    'config': '/home/di/code/openpipeline/src/query/cellxgene_census/config.vsh.yaml'
}
## VIASH END

def test_cellxgene_extract_metadata_expression(run_component, tmp_path):
    output_file = tmp_path / "output.h5ad"

    run_component([
        "--species", "homo_sapiens",
        "--obs_value_filter", "is_primary_data == True and cell_type_ontology_term_id in ['CL:0000136', 'CL:1000311', 'CL:0002616'] and suspension_type == 'cell'",
        "--output", output_file,
        "--obs_batch", "sex,sex",
        "--dataset_id", "test_dataset_id",
        "--dataset_name", "test_dataset_name",
        "--dataset_url", "https://test_dataset_url.com",
        "--dataset_reference", "test_dataset_reference",
        "--dataset_summary", "test_dataset_summary",
        "--dataset_description", "test_dataset_description",
        "--dataset_organism", "test_homo_sapiens",
    ])

    # check whether file exists
    assert os.path.exists(output_file), "Output file does not exist"

    adata = ad.read_h5ad(output_file)

    # check obs
    assert not adata.obs.empty, ".obs should not be empty"
    assert "is_primary_data" in adata.obs.columns
    assert np.all(adata.obs["is_primary_data"] == True)
    assert "cell_type_ontology_term_id" in adata.obs.columns
    assert "disease" in adata.obs.columns
    assert adata.n_obs > 10
    assert np.all([x in ["male+male", "female+female"] for x in adata.obs["batch"]])

    # check var
    assert "soma_joinid" in adata.var.columns
    assert "feature_id" in adata.var.columns

    # check uns
    assert adata.uns["dataset_id"] == "test_dataset_id", "Incorrect .uns['dataset_id']"
    assert adata.uns["dataset_name"] == "test_dataset_name", "Incorrect .uns['dataset_name']"
    assert adata.uns["dataset_url"] == "https://test_dataset_url.com", "Incorrect .uns['dataset_url']"
    assert adata.uns["dataset_reference"] == "test_dataset_reference", "Incorrect .uns['dataset_reference']"
    assert adata.uns["dataset_summary"] == "test_dataset_summary", "Incorrect .uns['dataset_summary']"
    assert adata.uns["dataset_description"] == "test_dataset_description", "Incorrect .uns['dataset_description']"
    assert adata.uns["dataset_organism"] == "test_homo_sapiens", "Incorrect .uns['dataset_organism']"


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
