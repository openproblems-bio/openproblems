from hashlib import md5
from typing import Callable, Any, Union, Optional

from io import BufferedIOBase
from os import PathLike
from urllib.error import HTTPError

from openproblems.datasets.callbacks._utils import _hash_function

PATH_CB_t = Callable[[Union[str, PathLike]], Any]
BUFF_CB_t = Callable[[BufferedIOBase], Any]
CB_t = Union[BUFF_CB_t, PATH_CB_t]
FP_t = Union[str, bytes, BufferedIOBase, PathLike]

from openproblems.datasets.exceptions import CallbackError  # TODO: move me
from openproblems.datasets.callbacks.exceptions import CallbackOSError, CallbackHTTPError  # TODO: or rather move me


class FunctionHasherMixin:

    def __init__(self):
        self._hash: Optional[int] = None

    def _create_hash(self, fn: Callable) -> int:
        if self._hash is None:
            hasher = md5()
            _hash_function(fn, hasher)
            self._hash = int(hasher.hexdigest(), base=16)

        return self._hash


class Callback(FunctionHasherMixin):

    def __init__(self, callback: CB_t):
        super().__init__()
        if not callable(callback):
            raise TypeError(type(callback))
        self._callback = callback

    def __hash__(self) -> int:
        return self._create_hash(self._callback)

    def __eq__(self, other) -> bool:
        return type(self) == type(other) and hash(self) == hash(other)

    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        try:
            return self._callback(*args, **kwargs)
        except OSError as e:
            raise CallbackOSError(self) from e
        except HTTPError as e:
            raise CallbackHTTPError(self) from e
        except Exception as e:
            raise CallbackError(self) from e

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}[callback={self._callback}]>"

    def __str__(self) -> str:
        return repr(self)
