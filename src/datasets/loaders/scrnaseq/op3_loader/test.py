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
    output_file = tmp_path / "output.h5ad"

    run_component([
        "--data_type", "sc",
        "--donor_id", "1",
        "--cell_type", "T cells",
        "--output", output_file,
        "--dataset_id", "op3_test",
        "--dataset_name", "OP3 Test Dataset",
        "--dataset_summary", "Test summary",
        "--dataset_description", "Test description",
    ])

    # check whether file exists
    assert os.path.exists(output_file), "Output file does not exist"

    adata = ad.read_h5ad(output_file)

    # check obs
    assert not adata.obs.empty, ".obs should not be empty"
    assert "donor_id" in adata.obs.columns
    assert "cell_type" in adata.obs.columns
    assert "perturbation" in adata.obs.columns
    assert adata.n_obs > 0

    # check var
    assert not adata.var.empty, ".var should not be empty"

    # check uns
    assert adata.uns["dataset_id"] == "op3_test", "Incorrect .uns['dataset_id']"
    assert adata.uns["dataset_name"] == "OP3 Test Dataset", "Incorrect .uns['dataset_name']"
    assert adata.uns["dataset_summary"] == "Test summary", "Incorrect .uns['dataset_summary']"
    assert adata.uns["dataset_description"] == "Test description", "Incorrect .uns['dataset_description']"
    assert adata.uns["dataset_organism"] == "Homo sapiens", "Incorrect .uns['dataset_organism']"

    # check layers
    assert "counts" in adata.layers, "counts layer should exist"


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))