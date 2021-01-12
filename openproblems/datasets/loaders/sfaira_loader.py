from typing import Optional, Any

from openproblems.datasets.callbacks.callback import Callback
from openproblems.datasets.loaders.loader import BaseLoader, URL_t


class SFAIRALoader(BaseLoader):
    def _download(self, url: URL_t, callback: Optional[Callback] = None, **kwargs) -> Any:
        raise NotImplementedError("SFAIRA loaders are not implemented.")

