from typing import Optional, Union, Sequence
from datetime import date
from dataclasses import dataclass, field
from openproblems.datasets.url import URL


@dataclass(frozen=True, order=False)
class Metadata:
    name: str
    author: str
    year: Union[int, date]  # TODO: warn if int?
    url: Union[str, URL, Sequence[Union[str, URL]]] = field(repr=False)
    description: Optional[str] = field(default=None, repr=False)

    def __post_init__(self):
        url = self.url
        if not isinstance(self.year, date):
            object.__setattr__(self, "url", date(year=self.year, month=1, day=1))

        if isinstance(url, str):
            url = (url,)
        url = tuple(u if isinstance(u, URL) else URL(u) for u in url)
        if len(url) == 1:  # no multiple files
            url = url[0]

        object.__setattr__(self, "url", url)
