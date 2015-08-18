"""Factory to load IAF files from an input StreamIAFFactory."""

from IAFFactory import IAFFactory


class StreamIAFFactory(IAFFactory):
    """Timeseries Factory for IAGA2002 formatted files loaded via a stream.
        normally either a single file, or stdio.

    Parameters
    ----------
    stream: file object
        io stream, normally either a file, or stdio

    See Also
    --------
    IAFFactory
    Timeseriesfactory
    """
    def __init__(self, stream, observatory=None, channels=None,
            type=None, interval=None, version1_flag=None):
        IAFFactory.__init__(self, None, observatory, channels,
            type, interval)
        self._stream = stream
        self._version1_flag = version1_flag

    def get_timeseries(self, starttime, endtime, observatory=None,
            channels=None, type=None, interval='minute'):
        """Implements get_timeseries

        Notes: Calls IAGFFactory.parse_string in place of
            IAFFactory.get_timeseries.
        """

        self._stream.seek(0)
        return IAFFactory.parse_string(self, self._stream.read(),
                self._version1_flag, interval)

    def put_timeseries(self, timeseries, starttime=None, endtime=None,
            channels=None, type=None, interval=None):
        """Implements put_timeseries
        """
        IAFFactory.write_file(self, self._stream, timeseries, channels)
