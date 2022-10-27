from .....tools.decorators import method
from .....tools.utils import check_version
from ..utils import aggregate_method_scores
import functools

import numpy as np
import pandas as pd

_random_method = functools.partial(
    method,
    paper_name="Random Events (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)


@_random_method(
    method_name="Random Events",
)
def random_events(adata, test=False, n_events=1000):
    adata.uns["ccc_pred"] = pd.DataFrame(
        {
            "ligand": np.random.choice(
                adata.uns["ligand_receptor_resource"]["ligand_genesymbol"], n_events
            ),
            "receptor": np.random.choice(
                adata.uns["ligand_receptor_resource"]["receptor_genesymbol"], n_events
            ),
            "source": np.random.choice(adata.obs["label"].cat.categories, n_events),
            "target": np.random.choice(adata.obs["label"].cat.categories, n_events),
            "score": np.random.uniform(0, 1, n_events),
        }
    )
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@_random_method(
    method_name="Random Events MAX",
)
def random_events_max(adata, test=False):
    adata = random_events(adata, test=test)
    adata.uns['ccc_pred'] = aggregate_method_scores(adata, how='max')

    return adata


@_random_method(
    method_name="Random Events SUM",
)
def random_events_sum(adata, test=False):
    adata = random_events(adata, test=test)
    adata.uns['ccc_pred'] = aggregate_method_scores(adata, how='sum')

    return adata


_true_method = functools.partial(
    method,
    paper_name="True Events (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)


@_true_method(
    method_name="True Events",
)
def true_events(adata, test=False):
    adata.uns["ccc_pred"] = adata.uns["ccc_target"].rename(
        {"response": "score"}, axis=1
    )
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata

@_true_method(
    method_name="True Events SUM",
)
def true_events_max(adata, test=False):
    adata = true_events(adata, test=test)
    adata.uns['ccc_pred'] = aggregate_method_scores(adata, how='max')

    return adata


@_random_method(
    method_name="True Events SUM",
)
def true_events_sum(adata, test=False):
    adata = true_events(adata, test=test)
    adata.uns['ccc_pred'] = aggregate_method_scores(adata, how='sum')

    return adata
