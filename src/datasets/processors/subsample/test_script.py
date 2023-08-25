import sys
import os
import pytest
import anndata as ad

## VIASH START
meta = {
    "resources_dir": "resources_test/common"
}
## VIASH END

input_path = f"{meta['resources_dir']}/pancreas/dataset.h5ad"
input = ad.read_h5ad(input_path)

def test_even_sampling(run_component):
    output_path = "output.h5ad"
    run_component([
        "--input", input_path,
        "--output", output_path,
        "--even",
        "--seed", "123",
        "--n_obs", "100",
        "--n_vars", "120"
    ])

    # Checking whether file exists
    assert os.path.exists(output_path), "Output file not found"

    # Check that test output fits expected API
    output = ad.read_h5ad(output_path)

    assert output.n_obs <= 100, "n_obs should be <= 100"
    assert output.n_vars <= 120, "n_vars should be <= 100"


def test_keep_functionality(run_component):
    output_path = "output.h5ad"

    # keep_features = list(input.var_names[:10])
    # use genes with high enough expression
    keep_features = ["ANP32E", "CBX5", "HMGB2"]

    run_component([
        "--input", input_path,
        "--keep_celltype_categories", "acinar:beta",
        "--keep_batch_categories", "celseq:inDrop4:smarter",
        "--keep_features", ":".join(keep_features),
        "--output", output_path,
        "--seed", "123"
    ])

    # Checking whether file exists
    assert os.path.exists(output_path), "Output file not found"

    # Check that test output fits expected API
    output = ad.read_h5ad(output_path)

    assert output.n_obs <= 500, "n_obs should be <= 500"
    assert output.n_vars <= 500, "n_vars should be <= 500"
    for feat in keep_features:
        assert feat in output.var_names, f"{feat} should be in output.var_names"

if __name__ == '__main__':
    sys.exit(pytest.main([__file__, "--capture=no"], plugins=["viashpy"]))
