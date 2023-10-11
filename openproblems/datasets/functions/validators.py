from openproblems.datasets.functions._utils import _validate_types, _get_types, _is_function
from openproblems.datasets.functions.base import FunctionBase
from anndata import AnnData
from abc import ABC
from typing import Tuple, Any, Type, Union, Iterable, Mapping, Optional, Callable, Sequence
from scipy.sparse import spmatrix
from pandas.core.dtypes.common import infer_dtype_from_object
from types import MappingProxyType
from enum import Enum
import pandas as pd
import numpy as np

Type_t = Tuple[Union[Type, str,], ...]

__all__ = ["ObjectTypeValidator", "ColumnTypeValidator", "IsAnnData"]


class AttributeLoc(str, Enum):
    OBS = "obs"
    OBSM = "obsm"
    VAR = "var"
    VARM = "varm"
    OBSP = "obsp"
    UNS = "uns"


class ValidatorFunction(FunctionBase, ABC):
    pass


class ObjectTypeValidator(ValidatorFunction):
    def __init__(self, types: Union[Type, Tuple[Type, ...]]):
        types = tuple(types) if isinstance(types, Iterable) else (types,)
        _validate_types(types)

        self._valid_types = types

    def __call__(self, obj: Any) -> Any:
        if not isinstance(obj, self._valid_types):
            raise TypeError(type(obj))
        return obj


class ColumnTypeValidator(ValidatorFunction):
    def __init__(self, types: Union[Sequence[Any], Mapping[Any, Optional[Type]]],
                 location: Union[str, AttributeLoc] = AttributeLoc.OBS,
                 exact: bool = True):
        location = AttributeLoc(location)

        if not isinstance(types, Mapping):
            types = {typ: None for typ in types}
        _validate_types(types.values(), allow_none=True, valid_types=(str, type, Callable))

        self._type_map = MappingProxyType({k: _get_types(v) for k, v in types.items()})
        self._location = location
        self._exact = exact

    def __call__(self, adata: AnnData) -> AnnData:
        haystack = getattr(adata, self._location)
        for colname, expected in self._type_map.items():
            if expected is None:
                continue

            needle = haystack[colname]

            if _is_function(expected):
                # warning: this is called on the object to support pandas checks
                res = expected(needle)
            else:
                if isinstance(needle, pd.Series):
                    actual = infer_dtype_from_object(needle.dtype)
                elif isinstance(needle, (np.ndarray, spmatrix)):
                    actual = infer_dtype_from_object(needle)
                else:
                    actual = type(needle)
                res = actual == expected if self._exact else issubclass(actual, expected)

            if not res:
                raise TypeError(colname, expected, actual)  # TODO: nicer message

        return adata


class IsAnnData(ObjectTypeValidator):
    def __init__(self):
        super().__init__(AnnData)
