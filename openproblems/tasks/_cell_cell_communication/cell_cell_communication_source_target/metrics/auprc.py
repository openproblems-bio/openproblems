from .....tools.decorators import metric
from ..._common.metrics import auprc as _auprc
from ..api import MERGE_KEYS


@metric(**_auprc.metadata)
def auprc(adata):
    return _auprc(adata, merge_keys=MERGE_KEYS)
