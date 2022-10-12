from .....tools.decorators import metric
from ..._common.metrics import prauc as _prauc
from ..api import MERGE_KEYS


@metric(**_prauc.metadata)
def prauc(adata):
    return _prauc(adata, merge_keys=MERGE_KEYS)
