"""Specific tests for the dimensionality_reduction task"""
import numpy as np
import openproblems
import openproblems.tasks._cell_cell_communication._common as common
import openproblems.tasks._cell_cell_communication._common.api
import openproblems.tasks._cell_cell_communication._common.utils
import os
import pandas as pd
import tempfile
import unittest
import utils.docker
import utils.git

TASK = openproblems.tasks.cell_cell_communication_source_target


class TestApi(unittest.TestCase):
    def test_assert_is_subset(self):
        assert (
            common.api.assert_is_subset(["a"], ["a", "b", "c"], prop_missing_allowed=0)
            is None
        )
        assert (
            common.api.assert_is_subset(
                ["a", "b", "c", "d"], ["a", "b", "c"], prop_missing_allowed=0.25
            )
            is None
        )
        self.assertRaisesRegex(
            AssertionError,
            r"test_subset is not a subset of test_superset\. d missing from"
            r" test_superset",
            common.api.assert_is_subset,
            ["a", "b", "c", "d"],
            ["a", "b", "c"],
            superset_name="test_superset",
            subset_name="test_subset",
            prop_missing_allowed=0,
        )
        self.assertRaisesRegex(
            AssertionError,
            r"Allowed proportion \(0.24\) of missing test_subset elements exceeded"
            r" \(0\.25\)\. d missing from test_superset",
            common.api.assert_is_subset,
            ["a", "b", "c", "d"],
            ["a", "b", "c"],
            superset_name="test_superset",
            subset_name="test_subset",
            prop_missing_allowed=0.24,
        )

    def test_map_gene_symbols(self):
        adata = common.api.sample_dataset(TASK.api.MERGE_KEYS)
        index = adata.var.index.to_numpy()
        index[0] = "many_to_one_1"
        index[1] = "many_to_one_2"
        index[2] = "one_to_many"
        index[3] = "one_to_one"
        index[4] = "one_to_none"
        adata.var.index = index
        map_df = pd.DataFrame(
            {
                "gene": [
                    "one_from_many",
                    "one_from_many",
                    "many_from_one_1",
                    "many_from_one_2",
                    "many_from_one_3",
                    "one_from_one",
                    "one_from_none",
                ],
                "alias": np.concatenate(
                    [adata.var.index[[0, 1, 2, 2, 2, 3]], ["none_to_one"]]
                ),
            }
        )
        with tempfile.TemporaryDirectory() as tempdir:
            map_filename = os.path.join(tempdir, "gene_map.csv")
            map_df.to_csv(map_filename)
            adata_mapped = common.utils.map_gene_symbols(adata, map_filename)
        self.assertEqual(adata_mapped.shape[0], adata.shape[0])
        self.assertEqual(adata_mapped.shape[1], adata.shape[1] + 1)
        adata_mapped.X = adata_mapped.X.tocsr()
        np.testing.assert_array_equal(
            adata[:, ["one_to_one"]].X.toarray(),
            adata_mapped[:, ["one_from_one"]].X.toarray(),
        )
        np.testing.assert_array_equal(
            adata[:, ["one_to_none"]].X.toarray(),
            adata_mapped[:, ["one_to_none"]].X.toarray(),
        )
        np.testing.assert_array_equal(
            adata[:, ["one_to_many", "one_to_many", "one_to_many"]].X.toarray(),
            adata_mapped[
                :, ["many_from_one_1", "many_from_one_2", "many_from_one_3"]
            ].X.toarray(),
        )
        np.testing.assert_array_equal(
            adata[:, ["many_to_one_1", "many_to_one_2"]].X.sum(axis=1).A,
            adata_mapped[:, ["one_from_many"]].X.toarray(),
        )
        self.assertNotIn("one_from_none", adata_mapped.var.index)


@utils.docker.docker_test(image=TASK.metrics.odds_ratio.metadata["image"])
def test_odds_ratio_no_match():  # pragma: nocover
    import numpy as np

    task = openproblems.tasks.cell_cell_communication_source_target
    metric = task.metrics.odds_ratio

    openproblems.log.debug(
        "Testing {} metric from {} task".format(metric.__name__, task.__name__)
    )

    adata = task.api.sample_dataset()

    # check expected output
    adata = task.api.sample_method(adata)
    m = metric(adata, top_prop=0.4)
    assert np.issubdtype("float64", m)
    assert m == 0.7

    # force perfect score
    adata = task.methods.true_events(adata)
    m = metric(adata, top_prop=0.4)
    assert m == 1

    # force exception
    adata.uns["ccc_target"]["response"] = 0
    m = metric(adata, top_prop=0.4)
    assert m is np.nan
