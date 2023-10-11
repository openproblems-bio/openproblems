from enum import Enum, auto, IntEnum


class NoValue(Enum):
    def __repr__(self):
        return f"<{self.__class__.__name__}.{self.name}>"


class Event(NoValue):
    CACHE_LOAD = auto()
    CACHE_STORE = auto()
    DOWNLOAD = auto()
    POST_DOWNLOAD = auto()
    PROCESS = auto()
    VERIFY = auto()
    SAMPLE = auto()

    # currently unused
    EVENT_GATHER = auto()
    EVENT_GATHER_TRANSFORM = auto()


class AutoPriority(IntEnum):
    def __new__(cls):
        value = len(cls.__members__) + 4 + 1  # 4 is just an offset so that we don't start from 0
        obj = int.__new__(cls, 2 ** value)
        obj._value_ = 2 ** value
        return obj


class EventPriority(AutoPriority):
    ULTRA = ()
    HIGH = ()
    MEDIUM = ()
    LOW = ()
