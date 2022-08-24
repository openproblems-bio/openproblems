from ....data.sample import load_sample_data

import numbers
import numpy as np
import pandas as pd
import scanpy as sc


def assert_is_subset(subset, superset, subset_name="subset", superset_name="superset"):
    """Assert `np.all(np.isin(subset, superset))` with a more readable error message"""
    is_missing = ~np.isin(subset, superset)
    if np.any(is_missing):
        x_missing = ",".join([x for x in subset[is_missing]])
        raise AssertionError(
            f"{subset_name} elements {x_missing} missing from {superset_name}"
        )


# Helper function to split complex subunits
def flatten_complex_subunits(entities):
    return np.unique(
        np.array([
            entity for entity_complex in entities for entity in entity_complex.split("_")
        ])
    )


def check_dataset(adata, merge_keys):
    """Check that dataset output fits expected API."""
    assert "label" in adata.obs
    assert "ccc_target" in adata.uns

    # check target organism
    assert "target_organism" in adata.uns
    assert isinstance(adata.uns["target_organism"], numbers.Integral)

    # check response columns
    assert "score" not in adata.uns["ccc_target"]
    assert "response" in adata.uns["ccc_target"]
    assert np.issubdtype(adata.uns["ccc_target"]["response"].dtype, int)
    assert np.all(np.isin(adata.uns["ccc_target"]["response"], [0, 1]))

    # check against resource
    if "ligand" in merge_keys or "receptor" in merge_keys:
        # verify resource is correct
        # TODO(@dbdimitrov): check this is right
        assert "ligand_receptor_resource" in adata.uns
        assert "ligand_receptor_resource" in adata.uns
        assert "receptor_genesymbol" in adata.uns["ligand_receptor_resource"]
        assert "ligand_genesymbol" in adata.uns["ligand_receptor_resource"]
        assert_is_subset(
            flatten_complex_subunits(
                adata.uns["ligand_receptor_resource"]["receptor_genesymbol"]
            ),
            adata.var.index,
            "resource receptor names",
            "gene names",
        )
        assert_is_subset(
            flatten_complex_subunits(
                adata.uns["ligand_receptor_resource"]["ligand_genesymbol"]
            ),
            adata.var.index,
            "resource ligand names",
            "gene names",
        )

    # check merge keys
    if "source" in merge_keys:
        assert "source" in adata.uns["ccc_target"]
        assert_is_subset(
            adata.uns["ccc_target"]["source"].unique(),
            adata.obs["label"].cat.categories,
            "source cell",
            "cell types",
        )
    if "target" in merge_keys:
        assert "target" in adata.uns["ccc_target"]
        assert_is_subset(
            adata.uns["ccc_target"]["target"].unique(),
            adata.obs["label"].cat.categories,
            "target cell",
            "cell types",
        )

    if "receptor" in merge_keys:
        # verify target receptors are in resource
        assert "receptor" in adata.uns["ccc_target"]
        assert_is_subset(
            adata.uns["ccc_target"]["receptor"].unique(),
            np.unique(adata.uns["ligand_receptor_resource"]["receptor_genesymbol"]),
            "target receptor names",
            "resource receptor names",
        )
    if "ligand" in merge_keys:
        # verify target ligands are in resource
        assert "ligand" in adata.uns["ccc_target"]
        assert_is_subset(
            adata.uns["ccc_target"]["ligand"].unique(),
            np.unique(adata.uns["ligand_receptor_resource"]["ligand_genesymbol"]),
            "target ligand names",
            "resource ligand names",
        )

    return True


def check_method(adata, merge_keys):
    """Check that method output fits expected API."""
    assert "ccc_pred" in adata.uns

    # check response columns
    assert "response" not in adata.uns["ccc_pred"]
    assert "score" in adata.uns["ccc_pred"]
    assert np.all(np.isreal(adata.uns["ccc_pred"]["score"]))

    # check merge keys
    if "ligand" in merge_keys:
        assert "ligand" in adata.uns["ccc_pred"]
        assert_is_subset(
            adata.uns["ccc_pred"]["ligand"].unique(),
            np.unique(adata.uns["ligand_receptor_resource"]["ligand_genesymbol"]),
            "predicted ligand names",
            "resource ligand names",
        )
    if "receptor" in merge_keys:
        assert "receptor" in adata.uns["ccc_pred"]
        assert_is_subset(
            adata.uns["ccc_pred"]["receptor"].unique(),
            np.unique(adata.uns["ligand_receptor_resource"]["receptor_genesymbol"]),
            "predicted receptor names",
            "resource receptor names",
        )
    if "source" in merge_keys:
        assert "source" in adata.uns["ccc_pred"]
        assert_is_subset(
            adata.uns["ccc_pred"]["source"].unique(),
            adata.obs["label"].cat.categories,
            "source cell",
            "cell types",
        )
    if "target" in merge_keys:
        assert "target" in adata.uns["ccc_pred"]
        assert_is_subset(
            adata.uns["ccc_pred"]["target"].unique(),
            adata.obs["label"].cat.categories,
            "target cell",
            "cell types",
        )

    return True


def sample_dataset(merge_keys):
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

    # generate target interactions
    adata.uns["ccc_target"] = pd.DataFrame(
        {
            "response": np.random.binomial(1, 0.2, 50),
            "ligand": np.random.choice(adata.var.index, 50),
            "receptor": np.random.choice(adata.var.index, 50),
            "source": np.random.choice(list(set(adata.obs.label)), 50),
            "target": np.random.choice(list(set(adata.obs.label)), 50),
        }
    )
    # subset columns
    adata.uns["ccc_target"] = adata.uns["ccc_target"][["response"] + merge_keys]

    # assign to human prior knowledge
    adata.uns["target_organism"] = 9606
    n_complexes = 5
    n_genes = len(adata.var.index)
    ligand_complexes = [
        "_".join(np.random.choice(adata.var.index, 2)) for _ in
        range(n_complexes)
    ]
    receptor_complexes = [
        "_".join(np.random.choice(adata.var.index, 2)) for _ in
        range(n_complexes)
    ]
    adata.uns["ligand_receptor_resource"] = pd.DataFrame(
        {
            "ligand_genesymbol": ligand_complexes + np.random.choice(
                adata.var.index, n_genes, replace=False).tolist(),
            "receptor_genesymbol": receptor_complexes + np.random.choice(
                adata.var.index, n_genes, replace=False).tolist(),
        }
    )

    return adata


def sample_method(adata, merge_keys):
    """Create sample method output for testing metrics in this task."""
    row_num = 500
    np.random.seed(1234)

    ligand_msk = ~adata.uns["ligand_receptor_resource"]["ligand_genesymbol"].isin(
        adata.var.index
    )
    receptor_msk = ~adata.uns["ligand_receptor_resource"]["receptor_genesymbol"].isin(
        adata.var.index
    )
    msk = np.logical_not(ligand_msk | receptor_msk)
    # keep only plausible interactions
    resource = adata.uns["ligand_receptor_resource"][msk]

    df = pd.DataFrame(np.random.random((row_num, 1)), columns=["score"])
    df["source"] = np.random.choice(np.unique(adata.obs[["label"]]), row_num)
    df["target"] = np.random.choice(np.unique(adata.obs[["label"]]), row_num)
    df["ligand"] = np.random.choice(
        np.unique(resource["ligand_genesymbol"].values), row_num
    )
    df["receptor"] = np.random.choice(
        np.unique(resource["receptor_genesymbol"].values), row_num
    )
    # subset columns
    df = df[["score"] + merge_keys]

    adata.uns["ccc_pred"] = df

    return adata
