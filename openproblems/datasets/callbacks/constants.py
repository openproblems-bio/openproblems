from enum import Enum, auto


class CallbackErrorPolicy(Enum):
    RAISE = auto()
    TRY_REST = auto()


class CallbackResolutionPolicy(Enum):
    PREFER_INFERRED = auto()
    PREFER_SUPPLIED = auto()
