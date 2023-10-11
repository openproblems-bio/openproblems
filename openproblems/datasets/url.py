from typing import Iterable, Tuple, Union, Optional, Callable


class URL(str):
    def __new__(cls, *args, fallback: Optional[Union[str, Iterable[str]]] = None, callback: Optional[Callable] = None):
        # TODO: type the callback properly
        res = super().__new__(cls, *args)
        if isinstance(fallback, (str, type(None))):
            fallback = (fallback,)
        # TODO: disallow URL in fallback?
        fallback = tuple(f for f in fallback if f != res)

        res._fallback = fallback
        res._callback = callback

        return res

    @property
    def callback(self) -> Optional[Callable]:
        return self._callback

    @property
    def fallback(self) -> Tuple[str, ...]:
        return self._fallback

    @property
    def next(self) -> "URL":
        if not len(self._fallback):
            raise RuntimeError("No more callbacks.")  # TODO: custom exc

        return URL(self.fallback[0], fallback=self.fallback[1:], callback=self.callback)
