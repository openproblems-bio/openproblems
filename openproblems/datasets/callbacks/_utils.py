import re
from functools import partial, singledispatchmethod
from hashlib import md5
from io import BufferedIOBase
from mimetypes import guess_type
from os import PathLike
from types import FunctionType, BuiltinFunctionType, MethodType, BuiltinMethodType
from typing import Callable, Union, Iterable, Any, Hashable, Optional


def _hash_function(fn: Callable, hasher: md5, *, encoding: str = "utf-8") -> None:
    def _is_pure_function(fn) -> bool:
        # `callable` not used because it would also catch our callback class + junk like classes
        return isinstance(fn, (FunctionType, BuiltinFunctionType, MethodType, BuiltinMethodType))

    def _hash_keys(ks: Union[str, Iterable[str]]) -> None:
        if isinstance(ks, str):
            ks = (ks,)

        for k in ks:
            hasher.update(k.encode(encoding))

    def _hash_values(vs: Iterable[Any]) -> None:
        for v in vs:
            if _is_pure_function(v):
                _hash_function(v, hasher, encoding=encoding)
            elif not isinstance(v, Hashable):
                hasher.update(repr(v).encode(encoding))
            else:
                h = hash(v)
                hasher.update(h.to_bytes((h.bit_length() + 7) // 8, byteorder="little"))

    try:
        code = fn.__code__
        hasher.update(code.co_code)
        _hash_keys(code.co_name)
        _hash_keys(code.co_filename)
        _hash_keys(code.co_names)
        _hash_values(code.co_consts)
    except AttributeError:
        if isinstance(fn, partial):
            # WARNING: stuff like partial(lambda x: 1, 3) and partial(lambda x: 1, x=3) will hash to different values
            _hash_function(fn.func, hasher, encoding=encoding)
            _hash_values(fn.args)
            _hash_keys(fn.keywords.keys())
            _hash_values(fn.keywords.values())
        else:
            hasher.update(repr(fn).encode(encoding))


class MimeType:
    try:
        from magic import Magic
        _m = Magic(mime=True, uncompress=True)
    except ImportError:
        _m = None

    _pat = re.compile(r"^(?P<type>\w+|\*)/(?P<tree>x(\.|-)|vnd\.|prs\.|.??)(?P<sub>\w+(-\w+)*|\*|)")

    # TODO: fix the type
    @classmethod
    def get(cls, obj: "FP_t") -> Optional[str]:
        mimetype = cls()._get(obj)
        if mimetype is None:
            return None

        match = cls._pat.match(mimetype)
        if match is None:
            # possibly malformed, don't care
            return mimetype

        # normalize
        return f"{match['type']}/{match['tree'].lower().replace('-', '.')}{match['sub']}"

    # TODO: fix the type
    @singledispatchmethod
    def _get(self, obj: "FP_t") -> Optional[str]:
        raise NotImplementedError(type(obj))

    @_get.register(str)
    @_get.register(PathLike)
    def _(self, path: PathLike) -> Optional[str]:
        try:
            if self._m is not None:
                return self._m.from_file(path)
        except OSError:
            pass
        except Exception:
            try:
                guess_type(path)
            except Exception:
                pass

    @classmethod
    @_get.register
    def _(self, fp: BufferedIOBase) -> Optional[str]:
        # if we were to surely use BytesIO, we could just get a part of the buffer and delegat to `bytes`
        try:
            # seek/tell will throw UnsupportedError, which is a subclass of OSError
            ix = fp.tell()
            fp.seek(0)
        except OSError:
            return

        try:
            if self._m is not None:
                return self._m.from_buffer(fp.read(8196))
        except Exception:
            pass
        finally:
            fp.seek(ix)

    @classmethod
    @_get.register
    def _(self, fp: bytes) -> Optional[str]:
        try:
            if self._m is not None:
                return self._m.from_buffer(fp)
        except Exception:
            pass