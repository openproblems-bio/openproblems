from ...data.sample import load_sample_data
from ...tools.decorators import dataset

import numpy as np
import pandas as pd
import scanpy as sc


def check_dataset(adata):
    """Check that dataset output fits expected API."""
    assert "label" in adata.obs
    assert "ccc_target" in adata.uns
    assert "response" in adata.uns["ccc_target"]
    assert np.issubdtype(adata.uns["ccc_target"]["response"].dtype, int)
    assert "target_organism" in adata.uns
    assert np.issubdtype(adata.uns["target_organism"], np.integer)
    return True


def check_method(adata):
    """Check that method output fits expected API."""
    assert "ccc_pred" in adata.uns
    # check if ligand-receptor and source-target (cell type) columns exist
    assert set(adata.uns["ccc_pred"]).issuperset(
        [
            "ligand",
            # 'receptor',
            "source",
            "target",
            "score",
        ]
    )

    return True


@dataset()
def sample_dataset():
    """Create a simple dataset to use for testing methods in this task."""
    adata = load_sample_data()

    # keep only the top 10 most variable
    sc.pp.highly_variable_genes(adata, n_top_genes=10)
    adata = adata[:, adata.var["highly_variable"]]
    # assign names to known interactions
    adata.var.index = [
        "LGALS9",
        "PTPRC",
        "LRP1",
        "CD47",
        "CD44",
        "COL1A1",
        "ADAM10",
        "SIRPA",
        "COL4A1",
        "THBS2",
    ]
    # transfer label
    adata.obs["label"] = adata.obs.cell_name

    adata.uns["ccc_target"] = pd.DataFrame(
        {
            "response": np.random.binomial(1, 0.2, 50),
            "ligand": np.random.choice(adata.var.index, 50),
            "target": np.random.choice(list(set(adata.obs.label)), 50),
        }
    )

    # assign to human prior knowledge
    adata.uns["target_organism"] = 9606
    return adata


def sample_method(adata):
    """Create sample method output for testing metrics in this task."""
    row_num = 10
    np.random.seed(1234)

    df = pd.DataFrame(np.random.random((row_num, 1)), columns=["score"])

    celltypes = list(map(pd.util.testing.rands, (3, 3, 4, 5)))
    lrs = list(map(pd.util.testing.rands, (4, 4, 4, 6)))

    df["source"] = np.random.choice(celltypes, row_num)
    df["target"] = np.random.choice(celltypes, row_num)
    df["ligand"] = np.random.choice(lrs, row_num)
    df["receptor"] = np.random.choice(lrs, row_num)

    adata.uns["ccc_pred"] = df
    return adata
