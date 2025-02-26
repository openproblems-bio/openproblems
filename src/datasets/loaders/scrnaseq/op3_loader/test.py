import sys
import os
import pytest
import anndata as ad
import numpy as np

## VIASH START
meta = {
    'resources_dir': './resources_test/',
    'executable': './target/docker/datasets/loaders/scrnaseq/op3_loader',
    'config': './src/datasets/loaders/scrnaseq/op3_loader/config.vsh.yaml'
}
## VIASH END

def test_op3_loader(run_component, tmp_path):
    """Test the OP3 loader."""
    output_file = tmp_path / "output.h5ad"

    run_component([
        "--donor_id", "1",
        "--cell_type", "T cells",
        "--min_cells", "3",
        "--min_genes", "200",
        "--dataset_id", "test_op3",
        "--dataset_name", "OP3 Test Dataset",
        "--dataset_summary", "Test summary for OP3 dataset",
        "--dataset_description", "Test description for OP3 dataset",
        "--output", output_file,
        "--output_compression", "gzip"
    ])

    # check whether file exists
    assert os.path.exists(output_file), "Output file does not exist"

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
    
    # check uns
    assert adata.uns["dataset_id"] == "test_op3", "Incorrect .uns['dataset_id']"
    assert adata.uns["dataset_name"] == "OP3 Test Dataset", "Incorrect .uns['dataset_name']"
    assert adata.uns["dataset_summary"] == "Test summary for OP3 dataset", "Incorrect .uns['dataset_summary']"
    assert adata.uns["dataset_description"] == "Test description for OP3 dataset", "Incorrect .uns['dataset_description']"


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))