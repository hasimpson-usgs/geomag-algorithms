
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
    """IAF writer.
    """

    def __init__(self, version = 2.2):
        self.headers = {
            'data_type': 'definitive'
        }
        self.data = numpy.array([], int)
        self.now = UTCDateTime()

    def write(self, out, timeseries, channels):
        """write timeseries to iaga file

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
        print 'start', starttime, endtime
        self._check_month(starttime, endtime)
        endtime = endtime + 60
        num_of_days = endtime.julday - stats['starttime'].julday
        print 'num_of_days', num_of_days
        channels = self._get_channels(channels, timeseries)
        for day in xrange(0, num_of_days):
            self._fill_headers(stats, channels,
                starttime + day * 86400, endtime)
            self._fill_data(timeseries, channels, starttime, endtime, day)
            # self._debug()
            self._write(out)

    def _get_channels(self, channels, timeseries):
        if channels is not None:
            return channels
        elif len(timeseries.select(channel='H')) > 0:
            return H_CHANNELS
        elif len(timeseries.select(channel='X')) > 0:
            return X_CHANNELS
        return None

    def _check_month(self, starttime, endtime):
        if endtime.month - starttime.month < 0:
            raise TimeseriesFactoryException(
                    'Time frame does not cover a full month "%s" "%s"' %
                            (starttime, endtime))

    def _write(self, out):
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

    def _fill_headers(self, stats, channels, starttime, endtime):
        """format headers for IAF file

        Parameters
        ----------
        stats: obspy.core.trace.stats
            holds the observatory metadata
        channels: array_like
            channels to be reported.

        Returns
        -------
        array_like
            an array containing formatted strings of header data.
        """

        """Parse header line.

        Adds value to ``self.headers``.
        """
        self.headers['station'] = stats['station'].rjust(4)
        print self.headers['station']
        print starttime
        self.headers['date'] = self._get_record_date(starttime)
        self.headers['geodetic_latitude'] = \
                int(float(stats['geodetic_latitude']) * 1000)
        self.headers['geodetic_longitude'] = \
                int(float(stats['geodetic_longitude']) * 1000)
        self.headers['elevation'] = int(stats['elevation'])
        self.headers['channels'] = stats['channels']
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

    def _get_record_date (self, starttime):
        return int('%04d%03d' % (starttime.year, starttime.day))

    def _get_publication_date(self):
        year = int(('%04d' % self.now.year)[2:4]) * 100
        month = self.now.month
        return str(year + month)

    def _get_version(self, version):
        return 3


    def _fill_data(self, timeseries, channels, starttime, endtime, day):
        """
        Parameters
        ----------
        timeseries : obspy.core.Stream
            stream containing traces with channel listed in channels
        channels : sequence
            list and order of channel values to output.
        """

        timeseries = self._get_slice(timeseries, starttime, day)
        stats = timeseries[0].stats
        self.data = numpy.array([], int)


        for channel, cnt in zip(channels, xrange(0,4)):
            print channel
            data = numpy.multiply(timeseries.select(
                    channel=channel)[0].data, 10)
            print 'length', len(data)
            if channel == 'G':
                if len(data) == \
                        len(data[numpy.isnan(data)]):
                    data[numpy.isnan(data)] = IAFParser.EIGHTS
            data[numpy.isnan(data)] = IAFParser.NINES
            self.data = numpy.append(self.data, data.astype(int))

        # We are assuming that we don't have hourly and daily values
        # for now, and calculating them. We could read this using an
        # input factory passed into IAFFactory if needed.
        data = numpy.empty(24)
        for channel, cnt in zip(channels, xrange(4,8)):
            hours = TimeseriesUtility.get_traces(
                    timeseries.select(channel=channel)[0])
            for hour in xrange(0,24):
                data[hour] = hours[hour].stats.statistics['average']
            data[:] = [x * 10 for x in data]
            self.data = numpy.append(self.data, data.astype(int))
        data = numpy.empty(1)
        for channel, cnt in zip(channels, xrange(8,12)):
            print 'about to call get_traces'
            days = TimeseriesUtility.get_traces(
                    timeseries.select(channel=channel)[0])
            data[0] = days[0].stats.statistics['average']
            self.data = numpy.append(self.data, data.astype(int))
        print self.data

        # Fix problem with slice not successfully dealing with a 3 hour span.
        print 'pre-trunc', timeseries.select(channel='K')[0].data
        data = timeseries.select(channel='K')[0].data[0:8].astype(int)
        data[numpy.isnan(data)] = IAFParser.NINES
        self.data = numpy.append(self.data, data)
        self.data = numpy.append(self.data, [0,0,0,0])
        print len(self.data)
        print self.data[5850:]

    def _get_slice(self, timeseries, starttime, day):
        start = starttime + day * 86400.0
        end = start + 86340.0
        print start, end
        return timeseries.slice(start, end)

    @classmethod
    def format(self, timeseries, channels):
        """Get an IAF formatted string.

        Calls write() with a StringIO, and returns the output.

        Parameters
        ----------
        timeseries : obspy.core.Stream
                Stream containing traces with channel listed in channels
        channels : sequence
                List and order of channel values to output.

        Returns
        -------
        unicode
          IAF formatted string.
        """
        out = StringIO()
        writer = IAFWriter()
        writer.write(out, timeseries, channels)
        return out.getvalue()
