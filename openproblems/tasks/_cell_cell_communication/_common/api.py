from ....data.sample import load_sample_data

import numbers
import numpy as np
import pandas as pd
import scanpy as sc

SAMPLE_RECEPTOR_NAMES = [
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


def assert_is_subset(
    subset,
    superset,
    subset_name="subset",
    superset_name="superset",
    prop_missing_allowed=0,
):
    """Assert `np.all(np.isin(subset, superset))` with a more readable error message"""
    subset = np.asarray(subset)
    is_missing = ~np.isin(subset, superset)
    prop_missing = np.sum(is_missing) / len(subset)

    if prop_missing > prop_missing_allowed:
        if prop_missing_allowed == 0:
            msg = f"{subset_name} is not a subset of {superset_name}. "
        else:
            msg = (
                f"Allowed proportion ({prop_missing_allowed}) of missing "
                f"{subset_name} elements exceeded ({prop_missing:.2f}). "
            )
        x_missing = ",".join([x for x in subset[is_missing]])
        raise AssertionError(msg + f"{x_missing} missing from {superset_name}")


# Helper function to split complex subunits
def flatten_complex_subunits(entities):
    return np.unique(
        np.array(
            [
                entity
                for entity_complex in entities
                for entity in entity_complex.split("_")
            ]
        )
    )


def check_dataset(adata, merge_keys):
    """Check that dataset output fits expected API."""
    assert "label" in adata.obs
    assert "ccc_target" in adata.uns

    assert "merge_keys" in adata.uns
    np.testing.assert_array_equal(adata.uns["merge_keys"], merge_keys)

    # check target organism
    assert "target_organism" in adata.uns
    assert isinstance(adata.uns["target_organism"], numbers.Integral)

    # check response columns
    assert "score" not in adata.uns["ccc_target"]
    assert "response" in adata.uns["ccc_target"]
    assert np.issubdtype(adata.uns["ccc_target"]["response"].dtype, int)
    assert np.all(np.isin(adata.uns["ccc_target"]["response"], [0, 1]))
    assert any(adata.uns["ccc_target"][merge_keys].duplicated()) is False

    # check against resource
    if "ligand" in merge_keys or "receptor" in merge_keys:
        assert "ligand_receptor_resource" in adata.uns
        assert "receptor_genesymbol" in adata.uns["ligand_receptor_resource"]
        assert "ligand_genesymbol" in adata.uns["ligand_receptor_resource"]
        assert "var_names_all" in adata.uns
        assert_is_subset(
            flatten_complex_subunits(
                adata.uns["ligand_receptor_resource"]["receptor_genesymbol"]
            ),
            adata.uns["var_names_all"],
            "resource receptor names",
            "gene names",
            0.1,
        )
        assert_is_subset(
            flatten_complex_subunits(
                adata.uns["ligand_receptor_resource"]["ligand_genesymbol"]
            ),
            adata.uns["var_names_all"],
            "resource ligand names",
            "gene names",
            0.1,
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


def check_method(adata, merge_keys, is_baseline=False):
    """Check that method output fits expected API."""
    assert "ccc_pred" in adata.uns

    # check response columns
    assert "response" not in adata.uns["ccc_pred"]
    assert "score" in adata.uns["ccc_pred"]
    assert np.all(np.isreal(adata.uns["ccc_pred"]["score"]))

    # Check if a single prediction is return for every merge_key combo
    assert (adata.uns["ccc_pred"].groupby(merge_keys).size() == 1).all()

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

    adata.uns["merge_keys"] = merge_keys

    # keep only the top 10 most variable
    sc.pp.highly_variable_genes(adata, n_top_genes=len(SAMPLE_RECEPTOR_NAMES))
    adata = adata[:, adata.var["highly_variable"]]
    # ensure we got the right number of genes
    adata = adata[:, : len(SAMPLE_RECEPTOR_NAMES)].copy()
    # hard-code var names to known interactions
    adata.var.index = adata.uns["var_names_all"] = SAMPLE_RECEPTOR_NAMES
    # transfer label
    adata.obs["label"] = adata.obs.cell_name

    # generate target interactions
    adata.uns["ccc_target"] = pd.DataFrame(
        {
            "ligand": np.random.choice(adata.var.index, 50),
            "receptor": np.random.choice(adata.var.index, 50),
            "source": np.random.choice(list(set(adata.obs.label)), 50),
            "target": np.random.choice(list(set(adata.obs.label)), 50),
        }
    )
    # drop duplicates
    adata.uns["ccc_target"] = adata.uns["ccc_target"].drop_duplicates(subset=merge_keys)
    # ensure positive response class is always present
    n_rows = adata.uns["ccc_target"].shape[0]
    response = np.zeros(n_rows, dtype=np.int64)
    response[0 : np.int(n_rows * 0.3)] = 1
    adata.uns["ccc_target"]["response"] = response
    # subset columns
    adata.uns["ccc_target"] = adata.uns["ccc_target"][["response"] + merge_keys]

    # assign to human prior knowledge
    adata.uns["target_organism"] = 9606
    n_complexes = 5
    n_genes = len(adata.var.index)
    ligand_complexes = [
        "_".join(np.random.choice(adata.var.index, 2)) for _ in range(n_complexes)
    ]
    receptor_complexes = [
        "_".join(np.random.choice(adata.var.index, 2)) for _ in range(n_complexes)
    ]
    adata.uns["ligand_receptor_resource"] = pd.DataFrame(
        {
            "ligand_genesymbol": np.concatenate(
                [
                    ligand_complexes,
                    np.random.choice(adata.var.index, n_genes, replace=False),
                ]
            ),
            "receptor_genesymbol": np.concatenate(
                [
                    receptor_complexes,
                    np.random.choice(adata.var.index, n_genes, replace=False),
                ]
            ),
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
