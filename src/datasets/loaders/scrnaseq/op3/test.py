
import os
import pytest
import anndata as ad
import numpy as np
import requests
import sys

## VIASH START
meta = {
    'resources_dir': './resources_test/',
    'executable': './target/docker/datasets/loaders/scrnaseq/op3_loader',
    'config': os.path.abspath('./src/datasets/loaders/scrnaseq/op3_loader/config.vsh.yaml')
}
## VIASH END

def test_op3_loader(run_component, tmp_path):
    """Test the OP3 loader."""
    input_file = tmp_path / "input.h5ad"  # Convert to string to be safe
    output_file = tmp_path / "output.h5ad"  # Convert to string to be safe

    # download file
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279945/suppl/GSE279945_sc_counts_processed.h5ad"
    response = requests.get(url)
    if response.status_code != 200:
        raise RuntimeError(f"Failed to download file from {url}, status code: {response.status_code}")

    # run loader
    run_component([
        "--input", str(input_file),
        "--var_feature_name", "index",
        "--donor_id", "1",
        "--cell_type", "T cells",
        "--dataset_id", "test_op3",
        "--dataset_name", "OP3 Test Dataset",
        "--dataset_summary", "Test summary for OP3 dataset",
        "--dataset_description", "Test description for OP3 dataset",
        "--dataset_url", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279945/suppl/GSE279945_sc_counts_processed.h5ad",
        "--output", str(output_file),
        "--output_compression", "gzip"
    ])

    # check whether file exists
    assert output_file.exists(), f"Output file {output_file} does not exist"

    adata = ad.read_h5ad(output_file)

    # check obs
    assert not adata.obs.empty, ".obs should not be empty"
    assert "donor_id" in adata.obs.columns
    assert np.all(adata.obs["donor_id"] == "1"), "Not all cells are from donor 1"
    assert "cell_type" in adata.obs.columns
    assert np.all(adata.obs["cell_type"] == "T cells"), "Not all cells are T cells"
    assert "perturbation" in adata.obs.columns
    assert adata.n_obs > 0

    # check var
    assert not adata.var.empty, ".var should not be empty"
    assert adata.n_vars > 0

    # check layers
    assert "counts" in adata.layers
    assert adata.X is not None, "X matrix should not be None"
    
    # check uns
    assert adata.uns["dataset_id"] == "test_op3", "Incorrect .uns['dataset_id']"
    assert adata.uns["dataset_name"] == "OP3 Test Dataset", "Incorrect .uns['dataset_name']"
    assert adata.uns["dataset_summary"] == "Test summary for OP3 dataset", "Incorrect .uns['dataset_summary']"
    assert adata.uns["dataset_description"] == "Test description for OP3 dataset", "Incorrect .uns['dataset_description']"

if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))