from .....tools.decorators import metric
from ..._common.metrics import odds_ratio as _odds_ratio
from ..api import MERGE_KEYS


@metric(**_odds_ratio.metadata)
def odds_ratio(adata):
    return _odds_ratio(adata, merge_keys=MERGE_KEYS)
