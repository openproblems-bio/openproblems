from urllib3.exceptions import HTTPError
from openproblems.datasets.exceptions import DatasetError


class CallbackOSError(DatasetError, OSError):
    pass


class CallbackHTTPError(DatasetError, HTTPError):
    pass
