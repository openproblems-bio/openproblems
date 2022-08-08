from ...data.sample import load_sample_data
from ...tools.decorators import dataset
from .utils import ligand_receptor_resource
from functools import reduce

import numbers
import numpy as np
import pandas as pd
import scanpy as sc

CCC_COLUMNS = {
    "ligand",
    "receptor",
    "source",
    "target",
}
CCC_TARGET_COLUMNS = CCC_COLUMNS.union({"response"})
CCC_PRED_COLUMNS = CCC_COLUMNS.union({"score"})


# Helper function to split complex subunits
def decomplexify(entities):
    entity_list = [np.chararray.split(x, sep="_").tolist() for x in entities]
    return reduce(lambda a, b: a + b, entity_list)


def check_dataset(adata):
    """Check that dataset output fits expected API."""
    assert "label" in adata.obs
    assert "ccc_target" in adata.uns
    # require no superfluous columns
    assert set(adata.uns["ccc_target"]).issubset(CCC_TARGET_COLUMNS)
    # require at least one non-response column to be defined
    assert len(set(adata.uns["ccc_target"]).intersection(CCC_TARGET_COLUMNS)) >= 2

    assert "response" in adata.uns["ccc_target"]
    assert np.issubdtype(adata.uns["ccc_target"]["response"].dtype, int)
    assert np.all(np.isin(adata.uns["ccc_target"]["response"], [0, 1]))
    if "source" in adata.uns["ccc_target"].columns:
        assert np.all(
            np.isin(
                adata.uns["ccc_target"]["source"].unique(),
                adata.obs["label"].cat.categories,
            )
        )
    if "target" in adata.uns["ccc_target"].columns:
        assert np.all(
            np.isin(
                adata.uns["ccc_target"]["target"].unique(),
                adata.obs["label"].cat.categories,
            )
        )

    resource = ligand_receptor_resource(adata.uns["target_organism"])
    if "receptor" in adata.uns["ccc_target"].columns:
        assert np.all(
            np.isin(
                adata.uns["ccc_target"]["receptor"].unique(),
                np.unique(resource["receptor_genesymbol"]),
            )
        )
    if "ligand" in adata.uns["ccc_target"].columns:
        assert np.all(
            np.isin(
                adata.uns["ccc_target"]["ligand"].unique(),
                np.unique(resource["ligand_genesymbol"]),
            )
        )

    assert "target_organism" in adata.uns
    assert isinstance(adata.uns["target_organism"], numbers.Integral)

    return True


def check_method(adata):
    """Check that method output fits expected API."""
    assert "ccc_pred" in adata.uns
    # check if ligand-receptor and source-target (cell type) columns exist
    assert set(adata.uns["ccc_pred"]).issuperset(CCC_PRED_COLUMNS)
    assert "response" not in adata.uns["ccc_pred"]

    assert np.all(
        np.isin(decomplexify(adata.uns["ccc_pred"]["ligand"].unique()), adata.var.index)
    )
    assert np.all(
        np.isin(
            decomplexify(adata.uns["ccc_pred"]["receptor"].unique()), adata.var.index
        )
    )

    resource = ligand_receptor_resource(adata.uns["target_organism"])
    assert np.all(
        np.isin(
            adata.uns["ccc_pred"]["ligand"].unique(),
            np.unique(resource["ligand_genesymbol"]),
        )
    )
    assert np.all(
        np.isin(
            adata.uns["ccc_pred"]["receptor"].unique(),
            np.unique(resource["receptor_genesymbol"]),
        )
    )

    assert np.all(
        np.isin(
            adata.uns["ccc_pred"]["source"].unique(),
            adata.obs["label"].cat.categories,
        )
    )
    assert np.all(
        np.isin(
            adata.uns["ccc_pred"]["target"].unique(),
            adata.obs["label"].cat.categories,
        )
    )

    assert np.all(np.isreal(adata.uns["ccc_pred"]["score"]))

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
    row_num = 500
    np.random.seed(1234)

    resource = ligand_receptor_resource(adata.uns["target_organism"])
    ligand_msk = ~resource["ligand_genesymbol"].isin(adata.var.index)
    receptor_msk = ~resource["receptor_genesymbol"].isin(adata.var.index)
    msk = np.logical_not(ligand_msk | receptor_msk)
    # keep only plausible interactions
    resource = resource[msk]

    df = pd.DataFrame(np.random.random((row_num, 1)), columns=["score"])
    df["source"] = np.random.choice(np.unique(adata.obs[["label"]]), row_num)
    df["target"] = np.random.choice(np.unique(adata.obs[["label"]]), row_num)
    df["ligand"] = np.random.choice(
        np.unique(resource["ligand_genesymbol"].values), row_num
    )
    df["receptor"] = np.random.choice(
        np.unique(resource["receptor_genesymbol"].values), row_num
    )

    adata.uns["ccc_pred"] = df

    return adata
