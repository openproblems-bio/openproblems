from abc import ABC
from anndata import AnnData
from openproblems.datasets.events.constants import Event, EventPriority
from openproblems.datasets.events.event_manager import PipeEventManager
from openproblems.datasets.exceptions import EmptyListenerError
from openproblems.datasets.loaders.loader import RemoteURLDownloader
from openproblems.exceptions import SCOPException
from openproblems.datasets.functions.processors import NoopProcessor
from openproblems.datasets.functions.samplers import SizedSampler
from openproblems.datasets.functions.validators import IsAnnData


def _(error: Exception, fn):
    def wrapper(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except Exception as e:
            raise error from e

    return wrapper


# TODO: add the metaclass, checking all fields and converting some fields
class DatasetABC(ABC):
    METADATA = None

    _DOWNLOADER = RemoteURLDownloader()
    _PROCESSOR = NoopProcessor()
    _SAMPLER = SizedSampler(size=500, cells=True) >> SizedSampler(500, cells=False)
    _VERIFIER = IsAnnData()

    _event_manager = PipeEventManager()

    def __init_subclass__(cls):
        cls.__initialize__()

    @classmethod
    def __initialize__(cls):
        pass
        # TODO: NYI
        # cls._event_manager.register(Event.CACHE_LOAD, ..., id="fetch-cache-data", priority=0)
        # cls._event_manager.register(Event.CACHE_STORE, ..., id="store-cache-data", priority=0)

    # TODO: different design:
    # make this a classmethod, returning some wrapper around anndata with 1 method (sample)
    # move the sampler to that method (descriptor most likely to downloader - to allow re-downloader)
    # alt:
    # register a file cache AND a memory cache (emit_first) for load (prioritize memcache), but emit for store
    # alt: use singleton objects
    def load(self, **kwargs) -> AnnData:
        url = self.METADATA.url
        # TODO: see __initialize__
        # data = _(SCOPException("Pre download error."))(self._event_manager.emit_first(Event.CACHE_LOAD, url))
        data = None

        if data is None:
            data = _(SCOPException("Download error."), self._event_manager.emit_first)(Event.DOWNLOAD, url, **kwargs)
            if not isinstance(url, str):  # multiple urls, need to combine them
                # TODO: maybe let _process handle it?
                data = _(SCOPException("Post download error."), self._event_manager.emit_first)(Event.POST_DOWNLOAD, data)
            data = self._process(data)

        self._verify(data)

        return data

    def sample(self, data: AnnData) -> AnnData:
        try:
            return self._verify(self._event_manager.emit_first(Event.SAMPLE, data))
        except EmptyListenerError:  # TODO: NYI
            raise  # TODO
        except SCOPException:
            raise
        except Exception as e:
            raise SCOPException("Sampling error.") from e

    def _process(self, data: AnnData) -> AnnData:
        try:
            data = self._event_manager.emit(Event.PROCESS, data)  # custom, class-specific processors
            # TODO: see __initialize__
            # self._event_manager.emit_first(Event.CACHE_STORE, data)
        except EmptyListenerError:  # TODO: NYI
            return data
        except SCOPException:
            raise
        except Exception as e:
            raise SCOPException("Processing error") from e

        return data

    def _verify(self, data: AnnData) -> AnnData:
        try:
            self._event_manager.emit(Event.VERIFY, data)
        except EmptyListenerError:  # TODO: NYI
            pass
        except SCOPException:
            raise
        except Exception as e:
            raise SCOPException("Unknown verification error.") from e

        return data


# TODO: metaclass to verify attributes are present and/or transform them (URL class + Function class)
class Dataset(DatasetABC):

    @classmethod
    def __initialize__(cls):
        super().__initialize__()
        cls._event_manager.register(Event.DOWNLOAD, cls._DOWNLOADER.download, id="default-download", priority=EventPriority.MEDIUM)
        cls._event_manager.register(Event.PROCESS, cls._PROCESSOR, id="default-processor", priority=EventPriority.MEDIUM)
        cls._event_manager.register(Event.VERIFY, cls._VERIFIER, id="default-verify", priority=EventPriority.MEDIUM)
        cls._event_manager.register(Event.SAMPLE, cls._SAMPLER, id="default-sample", priority=EventPriority.MEDIUM)

    # TODO: should DL be combined - not really
    # TODO: but URLs should have fallbacks (nonrec) to try (might need emit_first_noerror)

