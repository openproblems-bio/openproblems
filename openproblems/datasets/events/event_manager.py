from abc import ABC, abstractmethod
from functools import total_ordering
from typing import Callable, Any, Optional, List
from collections import defaultdict
from openproblems.datasets.callbacks.callback import Callback
import joblib as jl
from openproblems.datasets.events.constants import Event, EventPriority


@total_ordering
class EventHandler(Callback):

    def __init__(self, callback, id: Any, priority: int):
        super().__init__(callback)

        self._id = id
        self._priority = priority

    @property
    def id(self) -> Any:
        return self._id

    @property
    def priority(self) -> int:
        return self._priority

    def __lt__(self, other: "EventHandler") -> bool:
        if not isinstance(other, EventHandler):
            return NotImplemented
        return self.priority < other.priority

    def __eq__(self, other: object) -> bool:
        return super().__eq__(other) and self.priority == other.priority

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}[id={self.id}, priority={self.priority}, callback={self._callback}]>"

    def __str__(self) -> str:
        return repr(self)


class EventManagerABC(ABC):

    def __init__(self):
        self._handlers = defaultdict(list)

    def register(self, event: Event, callback: Callable[[Any], Any], id: Any, priority: int) -> "EventManagerABC":
        event = Event(event)
        self._handlers[event].append(EventHandler(callback=callback, id=id, priority=priority))

        return self

    def deregister(self, event: Event, id: Any) -> "EventManagerABC":
        event = Event(event)
        self._handlers[event] = [h for h in self._handlers[event] if h.id != id]

        return self

    def prioritize(self, event: Event, id: Any, new_priority: int) -> "EventManagerABC":
        return NotImplemented

    def handlers(self, event: Event) -> List[EventHandler]:
        return sorted(self._handlers[event])

    @abstractmethod
    def emit(self, event: Event, *args, **kwargs):
        pass

    def emit_first(self, event: Event, *args, **kwargs):
        for handler in self.handlers(event):
            res = handler(*args, **kwargs)
            if res is not None:
                return res


class PipeEventManager(EventManagerABC):
    def emit(self, event: Event, *args, **kwargs):
        handlers = self.handlers(event)
        if not len(args):
            raise ValueError()

        if not len(handlers):
            return None

        result = handlers[0](*args, **kwargs)
        args = args[1:]

        for handler in handlers[1:]:
            result = handler(result, *args, **kwargs)

        return result


# TODO: not used
class ScatterGatherEventManager(EventManagerABC):
    def __init__(self, n_jobs: Optional[int] = None, backend: str = "loky"):
        super().__init__()

        if n_jobs is not None and n_jobs == 0:
            raise ValueError()

        self._n_jobs = n_jobs
        self._backend = backend

        self.register(Event.EVENT_GATHER, lambda _: _, id="gather-identity", priority=EventPriority.MEDIUM)
        self.register(Event.EVENT_GATHER_TRANSFORM, lambda _: _, id="gather-transform-identity", priority=EventPriority.MEDIUM)

    def emit(self, event: Event, *args, **kwargs):
        res = jl.Parallel(n_jobs=self._n_jobs, backend=self._backend)(jl.delayed(h)(*args, **kwargs)
                                                                      for h in self.handlers(event))

        res = self.emit_first(Event.EVENT_GATHER, res)
        return self.emit_first(Event.EVENT_GATHER_TRANSFORM, res)
