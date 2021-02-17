import queue
import threading


class NonBlockingStreamReader:
    """Filestream that can be read without blocking.

    Parameters
    ----------
    stream: the stream to read from.
        Usually a process' stdout or stderr.
    """

    @staticmethod
    def populateQueue(stream, queue):
        """Collect lines from 'stream' and put them in 'queue'."""
        while True:
            line = stream.readline()
            if line:
                queue.put(line)
            else:
                raise RuntimeError("Stream ended without EOF.")

    def __init__(self, stream):
        """Initialize class."""
        self.stream = stream
        self.queue = queue.Queue()

        self.thread = threading.Thread(
            target=self.populateQueue, args=(self.stream, self.queue)
        )
        self.thread.daemon = True
        self.thread.start()  # start collecting lines from the stream

    def readline(self, timeout=None):
        """Read one line from stream if available.

        Parameters
        ----------
        timeout : int, optional (default: None)
            Seconds to wait for output. If None, block until response.
        """
        try:
            return self.queue.get(block=timeout is not None, timeout=timeout)
        except queue.Empty:
            return None

    def read(self, timeout=None):
        """Read all available lines from stream.

        Parameters
        ----------
        timeout : int, optional (default: None)
            Seconds to wait for output. If None, block until EOF.
        """
        output = b""
        while True:
            next_line = self.readline(timeout=timeout)
            if next_line is None:
                return output
            output += next_line + b"\n"
