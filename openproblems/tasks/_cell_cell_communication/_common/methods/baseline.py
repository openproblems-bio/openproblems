from .....tools.decorators import method
from .....tools.utils import check_version

import numpy as np
import pandas as pd


@method(
    method_name="Random Events",
    paper_name="Random Events (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
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


@method(
    method_name="True Events",
    paper_name="True Events (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def true_events(adata, test=False):
    adata.uns["ccc_pred"] = adata.uns["ccc_target"].rename(
        {"response": "score"}, axis=1
    )
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata
