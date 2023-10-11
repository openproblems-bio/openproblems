from abc import ABC, abstractmethod
from typing import Optional, Union, Any, Mapping, Sequence
from functools import singledispatchmethod

from openproblems.datasets.callbacks.constants import CallbackResolutionPolicy, CallbackErrorPolicy
from openproblems.datasets.url import URL
from requests import Request, Session
from requests.adapters import HTTPAdapter
from urllib3.util import Retry
from urllib3.exceptions import HTTPError
from io import BytesIO, BufferedIOBase
from tqdm.auto import tqdm

from os import PathLike

from openproblems.datasets.callbacks.callback import Callback
from openproblems.datasets.callbacks._utils import MimeType
from openproblems.datasets.callbacks.registry import CallBackRegistry, FILE_REG, BUFF_REG
from openproblems.datasets.exceptions import CallbackError
import joblib as jl

URL_t = Union[str, URL, PathLike]


class BaseLoader(ABC):
    _REGISTRY: CallBackRegistry = None
    _EVENT_MANAGER = None
    _FALLBACK_EXCEPTIONS = ()

    def __init__(self,
                 cb_resolution_policy: CallbackResolutionPolicy = CallbackResolutionPolicy.PREFER_SUPPLIED,
                 cb_error_policy: CallbackErrorPolicy = CallbackErrorPolicy.TRY_REST):
        self._cb_resolution_policy = cb_resolution_policy
        self._cb_error_policy = cb_error_policy

    @singledispatchmethod
    def download(self, url: Union[str, PathLike], **kwargs) -> Any:
        raise NotImplemented(type(url))

    @download.register(PathLike)
    @download.register(URL)
    @download.register(str)
    def _(self, url: URL_t, **kwargs) -> Any:
        callback = url.callback if isinstance(url, URL) else None
        if callback is not None:
            callback = Callback(callback)

        try:
            return self._download(url, callback=callback, **kwargs)
        except self._FALLBACK_EXCEPTIONS:
            if isinstance(url, URL):
                return self.download(url.next, **kwargs)
            raise
        except CallbackError as e:
            return self._resolve_callback_error(url, error=e, **kwargs)

    # TODO: maybe be less strict
    @download.register(frozenset)
    @download.register(set)
    @download.register(list)
    @download.register(tuple)
    def _(self, url: Sequence[URL_t], **kwargs) -> Any:
        with jl.Parallel(n_jobs=kwargs.pop("n_jobs", None)) as par:
            # warning: self.download can cause recursion
            return tuple(par(jl.delayed(self.download)(u, **kwargs) for u in url))

    def _resolve_callback_error(self, url: URL_t, *, error: CallbackError, **kwargs):
        if self._cb_error_policy == CallbackErrorPolicy.RAISE:
            raise error

        if self._cb_error_policy == CallbackErrorPolicy.TRY_REST:
            failed_cb, *_ = error.args
            if not isinstance(failed_cb, Callback):
                raise TypeError(f"Expected the failed callback to be of type `Callback`, "
                                f"found `{type(failed_cb).__name__}`.") from error

            for callback in self._REGISTRY.sort_closest_to(failed_cb):
                try:
                    # TODO
                    return self._download(url, callback=callback, **kwargs)
                except CallbackError:
                    pass
                except self._FALLBACK_EXCEPTIONS:
                    raise error from None  # most likely not recoverable
                except Exception:
                    pass

            raise error

        raise NotImplementedError(self._cb_error_policy)

    def _resolve_callback_preference(self, url: Union[str, PathLike, BufferedIOBase],
                                     callback: Optional[Callback] = None) -> Callback:
        try:
            raise  # check if we're currently handling a CallbackError
        except CallbackError:
            assert callback is not None, "Resolution callback is not set."
            return callback  # if so, override any preference since it doesn't matter
        except RuntimeError:
            pass  # not handling an any exception

        mime_type = MimeType.get(url)

        if callback is None:
            return self._REGISTRY[mime_type]

        if mime_type in self._REGISTRY:
            if self._cb_resolution_policy == CallbackResolutionPolicy.PREFER_SUPPLIED:
                return callback
            if self._cb_resolution_policy == CallbackResolutionPolicy.PREFER_INFERRED:
                return self._REGISTRY[mime_type]

        raise NotImplementedError(self._cb_resolution_policy)

    @abstractmethod
    def _download(self, url: URL_t, callback: Optional[Callback] = None, **kwargs) -> Any:
        pass


class LocalfileDownloader(BaseLoader):
    _REGISTRY = FILE_REG
    _FALLBACK_EXCEPTIONS = OSError

    def _download(self, url: URL_t, callback: Optional[Callback] = None, **kwargs) -> Any:
        return self._resolve_callback_preference(url, callback=callback)(url, **kwargs)


class RemoteURLDownloader(BaseLoader):
    _REGISTRY = BUFF_REG
    _FALLBACK_EXCEPTIONS = HTTPError

    def __init__(self, retries: int = 3, timeout: int = 600, chunk_size: int = 8192):
        # TODO: expose policies?
        super().__init__()
        self._session = Session()
        self._retries = retries
        self._timeout = timeout
        self._chunk_size = chunk_size

        if self._retries > 0:
            adapter = HTTPAdapter(
                max_retries=Retry(
                    total=self._retries,
                    redirect=5,
                    allowed_methods=["HEAD", "GET", "OPTIONS"],
                    status_forcelist=[413, 429, 500, 502, 503, 504],
                    backoff_factor=1,
                )
            )
            self._session.mount("ftp://", adapter)
            self._session.mount("http://", adapter)
            self._session.mount("https://", adapter)

    def _download(self, url: Union[str, PathLike], callback: Optional[Callback] = None,
                  params: Optional[Mapping[str, str]] = None, **kwargs) -> Any:
        handle = BytesIO()
        req = self._session.prepare_request(
            Request(
                "GET",
                url,
                params=params,
                headers={"User-agent": "openproblems-user"},
            )
        )

        with self._session.send(req, stream=True, timeout=self._timeout) as resp:
            resp.raise_for_status()
            total = resp.headers.get("content-length", None)

            with tqdm(
                    unit="B",
                    unit_scale=True,
                    miniters=1,
                    unit_divisor=1024,
                    total=total if total is None else int(total),
                    disable=kwargs.pop("disable", False),
            ) as t:
                for chunk in resp.iter_content(chunk_size=self._chunk_size):
                    t.update(len(chunk))
                    handle.write(chunk)

                handle.flush()
                handle.seek(0)

        return self._resolve_callback_preference(handle, callback=callback)(handle, ** kwargs)
