from abc import ABC, abstractmethod
from functools import partial
from types import MappingProxyType
from typing import Mapping, Mapping as TMapping, Iterable, Iterator

import pandas as pd
from anndata import read_h5ad as aread
from openproblems.datasets.callbacks.callback import CB_t, Callback
from scanpy import read as sread


class CallBackRegistry(Mapping, ABC):
    def __init__(self, data: TMapping[str, CB_t] = MappingProxyType({})):
        super().__init__()
        self._container = MappingProxyType({k:(Callback(v) if not isinstance(v, Callback) else v)
                                            for k, v in data.items()})
        for k, v in self._container.items():
            if not callable(v):  # check `callable` rather than our type (more future-proof)
                raise TypeError(f"Callback `{v}` registered for MIME `{k}` is not callable.")

    def sort_closest_to(self, v: Callback) -> Iterable[Callback]:
        # TODO: make something smarter
        return (cb for cb in self.values() if cb != v)

    def __getitem__(self, k: str) -> Callback:
        if k not in self._container:
            return self.__missing__(k)
        return self._container[k]

    def __contains__(self, k: str) -> bool:
        return k in self._container

    def __missing__(self, k: str) -> Callback:
        return self._default_callback

    def __len__(self) -> int:
        return len(self._container)

    def __iter__(self) -> Iterator["str"]:
        return iter(self._container)

    def __repr__(self) -> str:
        return repr(self._container)

    def __str__(self) -> str:
        return str(self._container)

    @property
    @abstractmethod
    def _default_callback(self):
        pass


class FileCallbackRegistry(CallBackRegistry):
    @property
    def _default_callback(self):
        return sread


class BufferCallbackRegistry(CallBackRegistry):
    @property
    def _default_callback(self):
        return aread


FILE_REG = FileCallbackRegistry({
    "text/comma-separated-values": pd.read_csv,
    "text/csv": pd.read_csv,
    "text/tab-separated-values": partial(pd.read_csv, sep="\t"),
    "text/tsv": partial(pd.read_csv, sep="\t"),
    "application/x-hdf5": sread,  # does not support file-like objects, only paths
})
BUFF_REG = BufferCallbackRegistry({
    "text/comma-separated-values": pd.read_csv,
    "text/csv": pd.read_csv,
    "text/tab-separated-values": partial(pd.read_csv, sep="\t"),
    "text/tsv": partial(pd.read_csv, sep="\t"),
    "application/x-hdf5": aread,  # does support file-like objects
})