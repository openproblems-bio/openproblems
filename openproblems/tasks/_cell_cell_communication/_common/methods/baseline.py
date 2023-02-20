from .....tools.decorators import baseline_method
from .....tools.utils import check_version

import numpy as np
import pandas as pd


@baseline_method(method_name="Random Events", method_summary="TODO")
def random_events(adata, test=False, n_events=1000):
    rng = np.random.default_rng(seed=1)

    ccc_pred = pd.DataFrame(
        {
            "ligand": rng.choice(
                adata.uns["ligand_receptor_resource"]["ligand_genesymbol"], n_events
            ),
            "receptor": np.random.choice(
                adata.uns["ligand_receptor_resource"]["receptor_genesymbol"], n_events
            ),
            "source": rng.choice(adata.obs["label"].cat.categories, n_events),
            "target": rng.choice(adata.obs["label"].cat.categories, n_events),
            "score": rng.uniform(0, 1, n_events),
        }
    )
    ccc_pred = ccc_pred.loc[~ccc_pred[adata.uns["merge_keys"]].duplicated()]

    adata.uns["ccc_pred"] = ccc_pred
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(method_name="True Events", method_summary="TODO")
def true_events(adata, test=False):
    adata.uns["ccc_pred"] = adata.uns["ccc_target"].rename(
        {"response": "score"}, axis=1
    )
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata
