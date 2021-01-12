from operator import __and__, __or__
from functools import reduce
from typing import Sequence, Callable, Optional, Union, Any
from openproblems.datasets.functions.base import FunctionBase
from anndata import AnnData
from abc import ABC, abstractmethod
from copy import copy

import numpy as np
import joblib as jl

from pandas.core.dtypes.common import is_categorical_dtype


def _optimize_index_sampler(fn: Callable, samplers: Sequence["SamplerFunction"]) -> Union["IndexSampler", Sequence[
    "SamplerFunction"]]:
    index_samplers = [s for s in samplers if isinstance(s, IndexSampler)]
    other_samplers = [s for s in samplers if not isinstance(s, IndexSampler)]

    assert len(index_samplers) + len(other_samplers) == len(samplers)

    if len(index_samplers) <= 1:
        return samplers

    cell_index = [s for s in index_samplers if s._cells]
    gene_index = [s for s in index_samplers if not s._cells]

    if len(cell_index) >= 2:
        dummy = copy(cell_index[0])
        dummy._calculate_index = lambda self, data: reduce(fn, [s._calculate_cell_index(data) for s in cell_index])
        new_cell_index = [dummy]
    else:
        new_cell_index = cell_index

    if len(gene_index) >= 2:
        dummy = copy(gene_index[0])
        dummy._calculate_index = lambda self, data: reduce(fn, [s._calculate_cell_index(data) for s in gene_index])
        new_gene_index = [dummy]
    else:
        new_gene_index = [dummy]

    return new_cell_index + new_gene_index + other_samplers


class SamplerFunction(FunctionBase, ABC):
    def __init__(self, cells: Optional[bool] = True):
        self._cells = cells

    def __and__(self, other: "SamplerFunction"):
        if not isinstance(other, SamplerFunction):
            return NotImplemented

        lsampler = self._samplers if isinstance(self, CombinedSampler) else [self]
        rsampler = other._samplers if isinstance(other, CombinedSampler) else [other]

        return AndCombiner(lsampler + rsampler)

    def __or__(self, other: "SamplerFunction"):
        if not isinstance(other, SamplerFunction):
            return NotImplemented

        lsampler = self._samplers if isinstance(self, CombinedSampler) else [self]
        rsampler = other._samplers if isinstance(other, CombinedSampler) else [other]

        return OrCombiner(lsampler + rsampler)

    @abstractmethod
    def _sample(self, data: AnnData) -> AnnData:
        pass

    def _pre_sample(self, data: AnnData, copy: bool = True):
        if self._cells or self._cells is None:
            return data.copy() if copy else data

        if data.is_view:
            data = data.copy()

        return data.T

    def _post_sample(self, data: AnnData):
        if self._cells is None:
            return data

        data = data if self._cells else data.T
        return data.copy() if data.is_view else data

    def __call__(self, data: AnnData, copy: bool = True):
        return self._post_sample(self._sample(self._pre_sample(data)))


class CombinedSampler(SamplerFunction, ABC):
    def __init__(self, samplers: Sequence[SamplerFunction], *, combiner: Callable, n_jobs: Optional[int] = None):
        super().__init__(cells=None)
        self._samplers = samplers
        self._combiner = combiner
        self._n_jobs = n_jobs

    def _sample(self, data: AnnData) -> AnnData:
        with jl.Parallel(n_jobs=self._n_jobs, prefer="processes") as par:
            subsets = par(jl.delayed(sampler)(data, copy=True) for sampler in self._samplers)

        obs_index = reduce(self._combiner, [subset.obs_name for subset in subsets], initial=data.obs_names)
        var_index = reduce(self._combiner, [subset.var_names for subset in subsets], initial=data.var_names)

        return data[obs_index, :][:, var_index]


class AndCombiner(CombinedSampler):
    def __init__(self, sampler: Sequence[SamplerFunction]):
        super().__init__(sampler, combiner=__and__)


class OrCombiner(CombinedSampler):
    def __init__(self, sampler: Sequence[SamplerFunction]):
        super().__init__(sampler, combiner=__or__)


class IndexSampler(SamplerFunction):
    def __init__(self, cells: bool = True):
        super().__init__(cells)

    def _sample(self, data: AnnData) -> AnnData:
        # TODO: ignore invalid option?
        data._inplace_subset_obs(self._calculate_index(data))
        return data

    @abstractmethod
    def _calculate_index(self, data: AnnData):
        pass


class ValueSampler(IndexSampler):
    def __init__(self, key: Any, values: Sequence[Any], *, cells: bool = True):
        super().__init__(cells)

        self._key = key
        self._values = values

    def _calculate_index(self, data: AnnData):
        return np.isin(data.obs[self._key], self._values)


class SizedSampler(IndexSampler):
    def __init__(self, size: int, *, cells: bool = True):
        super().__init__(cells)
        self._size = size

    def _calculate_index(self, data: AnnData):
        return np.arange(self._size)


class RandomSampler(SizedSampler):
    def __init__(self, size: int, *, seed: Optional[int] = None, cells: bool = True):
        super().__init__(size=size, cells=cells)
        self._rs = np.random.RandomState(seed=seed)

    def _calculate_index(self, data: AnnData):
        return self._rs.choice(data.n_obs, size=self._size, replace=False)


class GroupBalancedSampler(SamplerFunction):
    def __init__(self, keys: Union[str, Sequence[str]], *, sampler: SamplerFunction, cells: bool = True):
        super().__init__(cells=cells)
        self._keys = list(sorted(set(keys)))
        # TODO: cells are wrong
        self._sampler = sampler

    def _sample(self, data: AnnData) -> AnnData:
        obs = data.obs[self._keys]
        for col in obs.columns:
            if not is_categorical_dtype(obs[col]):
                raise TypeError(col)

        gb = obs.groupby(self._keys)
        subsets = []

        for group in gb.groups.keys():
            subset = data[(data.obs[self._keys] == group).all(1)]
            subset = self._sampler(subset)

            if subset.n_obs:
                subsets.append(subset)
            else:
                print("TODO: Empty subset")

        if not subsets:
            raise RuntimeError("All subsets are empty.")

        res = subsets[0].concatenate(*subsets[1:], batch_key=None)
        for col in data.obs.columns:
            if col in res.obs:
                res.obs[col] = res.obs[col].astype(data.obs[col].dtype)
        # TODO: remove -{batch_id} from index and keep batch index (but avoid overwrites)

        res.uns = data.uns
        res.varm = data.varm
        res.varp = data.varp

        return res
