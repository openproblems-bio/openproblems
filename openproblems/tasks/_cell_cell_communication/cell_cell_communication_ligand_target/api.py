from ....tools.decorators import dataset
from .._common import api

import functools

MERGE_KEYS = ["ligand", "target"]

check_dataset = functools.partial(api.check_dataset, merge_keys=MERGE_KEYS)
check_method = functools.partial(api.check_method, merge_keys=MERGE_KEYS)
sample_method = functools.partial(api.sample_method, merge_keys=MERGE_KEYS)


@dataset()
def sample_dataset():
    return api.sample_dataset(merge_keys=MERGE_KEYS)
