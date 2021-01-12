from typing import Callable, Any
from anndata import AnnData
from abc import ABC, abstractmethod


# TODO: rename to callbacks
class FunctionBase(ABC):
    def __rshift__(self, other: "FunctionBase") -> "FComposition":
        from openproblems.datasets.functions.composition import FComposition

        if not isinstance(other, FunctionBase):
            return NotImplemented

        left_fns = self._funcs if isinstance(self, FComposition) else [self]
        right_fns = self._funcs if isinstance(other, FComposition) else [other]

        return FComposition(left_fns + right_fns)

    @classmethod
    def from_callable(cls, fn: Callable[[AnnData], Any]) -> "FunctionBase":
        if isinstance(fn, cls):
            return fn
        if not callable(fn):
            raise TypeError()
        # TODO
        if not cls._validate_signature(fn):
            raise ValueError()

        return SimpleFunction(fn)

    @classmethod
    def _validate_signature(cls, fn: Callable[[AnnData], Any]) -> bool:
        # TODO: introduce expected signature
        pass

    @abstractmethod
    def __call__(self, *args, **kwargs) -> Any:
        pass


class SimpleFunction(FunctionBase):
    def __init__(self, fn: Callable[[AnnData], Any]):
        self._fn = fn

    def __call__(self, data: AnnData, **kwargs) -> Any:
        return self._fn(data, **kwargs)
