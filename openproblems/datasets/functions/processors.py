from abc import ABC

from typing import Optional, Union, Type, Sequence, Any

from pandas.core.dtypes.common import is_integer_dtype, is_bool_dtype, is_object_dtype, is_string_dtype, \
    is_categorical_dtype, is_signed_integer_dtype, is_unsigned_integer_dtype, is_float_dtype
from scipy.sparse import spmatrix, csr_matrix
from openproblems.datasets.functions.base import FunctionBase
from openproblems.datasets.functions._utils import _maybe_to_sparse_matrix, _maybe_convert_dtype
from anndata import AnnData
import scanpy as sc
import numpy as np


__all__ = ["NoopProcessor", "MakeUnique", "Sparsify", "UnifyDtypes", "Categorize"]


class ProcessorFunction(FunctionBase, ABC):
    pass


class NoopProcessor(ProcessorFunction):
    def __call__(self, adata: AnnData) -> AnnData:
        return adata


class MakeUnique(ProcessorFunction):
    def __call__(self, adata: AnnData) -> AnnData:
        adata.obs_names_make_unique()
        adata.var_names_make_unique()
        return adata


class RemoveEmpty(ProcessorFunction):
    def __call__(self, data: AnnData) -> AnnData:
        sc.pp.filter_cells(data)
        sc.pp.filter_genes(data)
        return data


class Sparsify(ProcessorFunction):
    def __init__(self, X: bool = True, layers: Optional[Sequence[str]] = None, raw: bool = True,
                 cls: Type[spmatrix] = csr_matrix):
        super().__init__()
        self._X = X
        self._raw = raw
        self._layers = layers
        self._cls = cls

    def __call__(self, data: AnnData) -> AnnData:
        if self._X:
            data.X = _maybe_to_sparse_matrix(data.X, cls=self._cls)

        if self._raw:
            if data.raw is not None:
                data.raw._X = _maybe_to_sparse_matrix(data.raw.X, cls=self._cls)
            else:
                pass  # TODO: log

        layers = data.layers.keys() if self._layers is None else self._layers

        for layer in layers:
            if layer in data.layers:
                data.layers[layer] = _maybe_to_sparse_matrix(data.layers[layer], cls=self._cls)
                pass
            else:
                pass  # TODO: log

        return data


class UnifyDtypes(ProcessorFunction):
    def __init__(self, signed: np.dtype = np.int32, unsigned: np.dtype = np.uint32, floating: np.dtype = np.float64):
        if not is_signed_integer_dtype(signed):
            raise TypeError()
        if not is_unsigned_integer_dtype(unsigned):
            raise TypeError()
        if not is_float_dtype(floating):
            raise TypeError()

        super().__init__()

        self._dtype_map = {np.signedinteger: signed, np.unsignedinteger: unsigned, np.floating: floating}

    def __call__(self, data: AnnData) -> AnnData:
        data.X = _maybe_convert_dtype(data.X, self._dtype_map)
        if data.raw is not None:
            data.raw._X = _maybe_convert_dtype(data.raw.X, self._dtype_map)

        for layer in data.layers:
            data[layer] = _maybe_convert_dtype(data.layers[layer], self._dtype_map)

        return data


class Categorize(ProcessorFunction):
    def __init__(self, key: Union[Any, Sequence[Any]]):
        if isinstance(key, str):
            key = (key,)
        if not isinstance(key, Sequence):
            raise TypeError()

        self._keys = sorted(key)

    def __call__(self, data: AnnData) -> AnnData:
        for key in self._keys:
            if key not in data.obs:
                continue  # TODO: log
            col = data.obs[key]

            if not is_categorical_dtype(col):
                if is_integer_dtype(col) or is_bool_dtype(col) or is_string_dtype(col) or is_object_dtype(col):
                    data.obs[key] = col.astype("category").values
                else:
                    continue  # TODO: log

            data.obs[key].cat.remove_empty_categories(inplace=True)

        return data
