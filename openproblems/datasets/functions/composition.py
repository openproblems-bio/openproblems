from functools import singledispatchmethod
from typing import Sequence, Any

from anndata import AnnData
from openproblems.datasets.functions.base import FunctionBase
from openproblems.datasets.functions.samplers import SamplerFunction
from openproblems.datasets.functions.processors import ProcessorFunction
from openproblems.datasets.functions.validators import ValidatorFunction


class FComposition(FunctionBase):
    def __init__(self, funcs: Sequence[FunctionBase]):
        # TODO: verify signatures
        assert len(funcs) >= 2
        self._funcs = funcs

    def __call__(self, data: AnnData, **kwargs) -> Any:
        for fn in self._funcs:
            data = self._handle_result(fn, data, fn(data, **kwargs))
        return data

    @singledispatchmethod
    def _handle_result(self, mixin: "FunctionBase", data: AnnData, result: Any) -> Any:
        raise NotImplementedError(mixin)

    @_handle_result.register(ProcessorFunction)
    @_handle_result.register(SamplerFunction)
    def _(self, _mixin: "FunctionBase", _data: AnnData, result: AnnData) -> Any:
        return result

    @_handle_result.register(ValidatorFunction)
    def _(self, mixin: "FunctionBase", data: AnnData, result: bool) -> Any:
        if not result:
            raise RuntimeError("Validator failed:", mixin)  # TODO: custom exception
        return data