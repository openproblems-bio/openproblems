from types import FunctionType, MethodType, BuiltinFunctionType, BuiltinMethodType
from typing import Iterable, Optional, Type, Union, Tuple, Callable, Mapping, Any

import builtins
import numpy as np
from pandas.core.dtypes.common import infer_dtype_from_object
from scipy.sparse import spmatrix, csr_matrix, issparse
from inspect import isclass


builtin_types = tuple(getattr(builtins, d) for d in dir(builtins) if isinstance(getattr(builtins, d), type))


def _is_function(fn):
    return isinstance(fn, FunctionType) or isinstance(fn, MethodType) or \
           isinstance(fn, BuiltinFunctionType) or isinstance(fn, BuiltinMethodType)


def _validate_types(types: Iterable[Optional[Type]], *,
                    allow_none: bool = False,
                    valid_types: Union[Type, Tuple[Type, ...]] = type) -> None:
    for typ in types:
        if typ is None:
            if not allow_none:
                raise ValueError()
            continue
        if not isinstance(typ, valid_types):
            raise TypeError(type(typ))


def _get_types(typ: Optional[Union[str, Callable, type]]) -> Union[Callable[[Any], bool], np.dtype, type]:
    if typ is None:
        return lambda _: True
    if _is_function(typ):
        return typ  # TODO and a warning: the callbacks are called on the object, not the dtype
    if typ in builtin_types:
        return typ
    if isclass(typ) and issubclass(typ, type):
        return typ
    return infer_dtype_from_object(typ)  # get a numpy-like dtype


def _maybe_to_sparse_matrix(arr: Union[np.ndarray, spmatrix], cls: Type[spmatrix] = csr_matrix) -> Union[np.ndarray, spmatrix]:
    if issparse(arr):
        if not isinstance(arr, cls):
            return cls(arr)

    if np.issubdtype(arr.dtype, np.integer):
        return cls(arr)

    # TODO: logg

    return arr


def _maybe_convert_dtype(arr: np.ndarray, dtype_map: Mapping[np.dtype, np.dtype]) -> np.ndarray:
    for old_type, new_type in dtype_map.items():
        if np.issubdtype(arr.dtype, old_type):
            return arr.astype(new_type, copy=False)

    return arr
