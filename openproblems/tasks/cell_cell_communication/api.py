from ...data.sample import load_sample_data
from ...tools.decorators import dataset
from random import choices
from random import seed

import numpy as np
import pandas as pd
import scanpy as sc


def check_dataset(adata):
    """Check that dataset output fits expected API."""
    assert "label" in adata.obs
    assert "bench" in adata.uns
    assert "response" in adata.uns["bench"]
    assert "target_organism" in adata.uns
    assert np.isreal(adata.uns["target_organism"])
    return True


def check_method(adata):
    """Check that method output fits expected API."""
    assert "ccc" in adata.uns
    # check if ligand-receptor and source-target (cell type) columns exist
    assert set(adata.uns["ccc"]).issuperset(
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

    adata.uns["bench"] = pd.DataFrame(
        {"response": np.random.binomial(1, 0.2, 50),
         "ligand": choices(adata.var.index, k=50),
         "target": choices(list(set(adata.obs.label)), k=50)}
    )

    # assign to human prior knowledge
    adata.uns["target_organism"] = 9606
    return adata


def sample_method(adata):
    """Create sample method output for testing metrics in this task."""
    row_num = 10
    seed(1234)

    df = pd.DataFrame(np.random.random((row_num, 1)), columns=["score"])

    celltypes = list(map(pd.util.testing.rands, (3, 3, 4, 5)))
    lrs = list(map(pd.util.testing.rands, (4, 4, 4, 6)))

    df["source"] = choices(celltypes, k=row_num)
    df["target"] = choices(celltypes, k=row_num)
    df["ligand"] = choices(lrs, k=row_num)
    df["receptor"] = choices(lrs, k=row_num)

    adata.uns["ccc"] = df
    return adata
