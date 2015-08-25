
from cStringIO import StringIO
import numpy
import struct
from ..TimeseriesFactoryException import TimeseriesFactoryException
from .. import TimeseriesUtility
from obspy.core import UTCDateTime
import IAFParser

X_CHANNELS = ['X', 'Y', 'Z', 'G']
H_CHANNELS = ['H', 'D', 'Z', 'G']


class IAFWriter(object):
    """IAFV283 writer.

    Attributes
    ----------
    headers : dict
        parsed IMFV283 headers.
    data : array_like
        all data values for a single day
    now : obspy.core.UTCDateTime
        current time
    version : int
        integer representing version
        0 - version 1.0
        1 - version 1.1
        2 - version 2.0
        3 - version 2.1
    """

    def __init__(self, version = 3):
        self.headers = {
            'data_type': 'definitive'
        }
        self.data = numpy.array([], int)
        self.now = UTCDateTime()
        self.version = version

    def write(self, out, timeseries, channels):
        """write timeseries to IMFV283 file

        Parameters
        ----------
        out: file object
            file object to be written to. could be stdout
        timeseries: obspy.core.stream
            timeseries object with data to be written
        channels: array_like
            channels to be written from timeseries object
        """
        stats = timeseries[0].stats
        starttime = stats['starttime']
        endtime = stats['endtime']
        self._check_month(starttime, endtime)
        endtime = endtime + 60
        num_of_days = endtime.julday - stats['starttime'].julday
        channels = self._get_channels(channels, timeseries)
        # each record is a single day
        for day in xrange(0, num_of_days):
            self._fill_headers(stats, channels, starttime + day * 86400)
            self._fill_data(timeseries, channels, starttime, day)
            self._write(out)

    def _check_month(self, starttime, endtime):
        """Make certain time frame is a month

        Parameters
        ----------
        starttime : obspy.core.UTCDateTime
            starttime of data to be written
        endtime : obspy.core.UTCDateTime
            endtime of data to be written

        Raises
        ------
        TimeseriesFactoryException
            if time frame is not exactly one month
        """
        # TODO might want to refactor this for multiple months
        end = endtime + 60
        if ((end.month - starttime.month ) != 1) or (end.day != 1):
            raise TimeseriesFactoryException(
                    'Time frame does not cover a full month "%s" "%s"' %
                            (starttime, endtime))

    def _get_channels(self, channels, timeseries):
        """get channels
        Parameters
        ----------
        channels : array_like
            array of channels passed in from controller
        timeseries : obspy.core.Stream
            timeseries data
        Returns
        -------
        array_like
            an array of channels
        """
        if channels is not None:
            return channels
        elif len(timeseries.select(channel='H')) > 0:
            return H_CHANNELS
        elif len(timeseries.select(channel='X')) > 0:
            return X_CHANNELS
        return None

    def _get_channels_headers(self, channels, stats):
        """get channels headers

        Parameters
        ----------
        channels : array_like
            array of channels
        stats : obspy.core.stats
            metadata

        Returns
        -------
        str
            a str to fill the channels field of the IAFV283 file
        """
        if 'channels' in stats:
            return stats['channels']
        else:
            if channels is not None:
                return ''.join(channels)
            else:
                raise TimeseriesFactoryException('No channels found')

    def _get_publication_date(self):
        """get publication date

        Returns
        -------
        str
            The current date as a string formated as the year and month YYmmm
        """
        year = int(('%04d' % self.now.year)[2:4]) * 100
        month = self.now.month
        return str(year + month)

    def _get_record_date (self, starttime):
        """get record data

        Parameters
        ----------
        starttime : obspy.core.UTCDateTime
            starttime of data

        Returns
        -------
        int
            The year and day of year as an int in the format YYYYddd
        """
        return int('%04d%03d' % (starttime.year, starttime.day))

    def _get_version(self, version):
        """get version

        Parameters
        ----------
        version : int
            version of the file
        Returns
        -------
            version of the file
        """
        return version or self.version

    def _fill_data(self, timeseries, channels, starttime, day):
        """
        Parameters
        ----------
        timeseries : obspy.core.Stream
            stream containing traces with channel listed in channels
        channels : sequence
            list and order of channel values to output.
        starttime : obspy.core.UTCDateTime
            starttime of the data
        day : int
            the day in the data we are working on
        Notes:
        Data consists of 13 channels/intervals. 4 minute channels,
            4 hour channels, 4 day channels, and 1 tri-hourly k channel
        """

        timeseries = self._get_slice(timeseries, starttime, day)
        stats = timeseries[0].stats
        self.data = numpy.array([], int)

        for channel, cnt in zip(channels, xrange(0,4)):
            if timeseries.select(channel=channel).count() == 0:
                raise TimeseriesFactoryException(
                    'Channel %s is required', )
            data = numpy.multiply(timeseries.select(
                    channel=channel)[0].data, 10)
            if channel == 'G':
                if len(data) == \
                        len(data[numpy.isnan(data)]):
                    data[numpy.isnan(data)] = IAFParser.EIGHTS
            data[numpy.isnan(data)] = IAFParser.NINES
            self.data = numpy.append(self.data, data.astype(int))

        # We are assuming that we don't have hourly and daily values
        # for now, and calculating them. We could read this using an
        # input factory passed into IAFFactory if needed.
        # f/G is NOT written out. Fill with
        data = numpy.empty(24)
        for channel, cnt in zip(channels, xrange(4,8)):
            hours = TimeseriesUtility.get_traces(
                    timeseries.select(channel=channel)[0])
            if channel == 'G':
                for hour in xrange(0,24):
                    data[hour] = IAFParser.NINES
            else:
                for hour in xrange(0,24):
                    data[hour] = hours[hour].stats.statistics['average']
                data[:] = [x * 10 for x in data]
            self.data = numpy.append(self.data, data.astype(int))
        data = numpy.empty(1)
        for channel, cnt in zip(channels, xrange(8,12)):
            days = TimeseriesUtility.get_traces(
                    timeseries.select(channel=channel)[0])
            data[0] = days[0].stats.statistics['average']
            self.data = numpy.append(self.data, data.astype(int))

        # Fix problem with slice not successfully dealing with a 3 hour span.
        data = timeseries.select(channel='K')[0].data[0:8].astype(int)
        data[numpy.isnan(data)] = IAFParser.NINES
        self.data = numpy.append(self.data, data)
        self.data = numpy.append(self.data, [0,0,0,0])

    def _fill_headers(self, stats, channels=None, starttime=None):
        """format headers for IAFV283 file

        Parameters
        ----------
        stats : obspy.core.trace.stats
            holds the observatory metadata
        channels : array_like
            channels to be reported.
        starttime : obspy.core.UTCDateTime
            starttime of the data
        """
        self.headers['station'] = stats['station'].rjust(4)
        self.headers['date'] = self._get_record_date(starttime)
        self.headers['geodetic_latitude'] = \
                int(float(stats['geodetic_latitude']) * 1000)
        self.headers['geodetic_longitude'] = \
                int(float(stats['geodetic_longitude']) * 1000)
        self.headers['elevation'] = int(stats['elevation'])
        self.headers['channels'] = self._get_channels_headers(
                channels, stats)
        self.headers['agency_name'] = stats['agency_name']
        self.headers['d_conversion'] = stats['d_conversion']
        self.headers['data_quality'] = stats['data_quality']
        self.headers['instrumentation'] = stats['instrumentation']
        self.headers['k_9'] = stats['k_9']
        self.headers['sensor_sampling_rate'] = \
                int(float(stats['sensor_sampling_rate']) * 1000.0)
        self.headers['sensor_orientation'] = stats['sensor_orientation']
        self.headers['publication_date'] = self._get_publication_date()
        self.headers['version'] = self._get_version(stats['version'])
        self.headers['empty'] = 0

    def _get_slice(self, timeseries, starttime, day):
        """get slice
        timeseries : obspy.core.stream
            timeseries data
        starttime : obspy.core.UTCDateTime
            starttime of the data
        day : int
            number of days from starttime
        """
        start = starttime + day * 86400.0
        end = start + 86340.0
        return timeseries.slice(start, end)

    def _write(self, out):
        """ write record to file

        Parameters
        ----------
        out : output pipe
            pre-opened write output pipe
        """
        record = struct.pack(IAFParser.IAFSTRING,
            self.headers['station'],
            self.headers['date'],
            self.headers['geodetic_latitude'],
            self.headers['geodetic_longitude'],
            self.headers['elevation'],
            self.headers['channels'],
            self.headers['agency_name'],
            self.headers['d_conversion'],
            self.headers['data_quality'],
            self.headers['instrumentation'],
            self.headers['k_9'],
            self.headers['sensor_sampling_rate'],
            self.headers['sensor_orientation'],
            self.headers['publication_date'],
            self.headers['version'],
            0,
            *self.data)
        out.write(record)
