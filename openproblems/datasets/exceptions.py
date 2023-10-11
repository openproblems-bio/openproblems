from openproblems.exceptions import SCOPException


class DatasetError(SCOPException):
    pass


# TODO: move me
class CallbackError(DatasetError, RuntimeError):
    pass


# TODO: move me
class EmptyListenerError(SCOPException, ValueError):
    pass
