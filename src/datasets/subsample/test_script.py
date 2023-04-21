import sys
import os
import pytest
import anndata as ad
import numpy as np

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

    print(">> Checking whether file exists")
    assert os.path.exists(output_path)

    print(">> Check that test output fits expected API")
    output = ad.read_h5ad(output_path)

    assert output.n_obs <= 100
    assert output.n_vars <= 120

    print(f"output: {output}", flush=True)


def test_keep_functionality(run_component):
    output_path = "output.h5ad"

    keep_features = list(input.var_names[:10])
    run_component([
        "--input", input_path,
        "--keep_celltype_categories", "acinar:beta",
        "--keep_batch_categories", "celseq:inDrop4:smarter",
        "--keep_features", ":".join(keep_features),
        "--output", output_path,
        "--seed", "123"
    ])

    print(">> Checking whether file exists")
    assert os.path.exists(output_path)

    print(">> Check that test output fits expected API")
    output = ad.read_h5ad(output_path)

    assert output.n_obs <= 500
    assert output.n_vars <= 500
    assert np.all([ f in output.var_names for f in keep_features])

    print(f"output: {output}", flush=True)

if __name__ == '__main__':
    sys.exit(pytest.main([__file__, "--capture=no"], plugins=["viashpy"]))










