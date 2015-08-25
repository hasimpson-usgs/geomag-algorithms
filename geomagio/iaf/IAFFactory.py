"""Factory that loads IAF Files."""

import obspy.core
import os
import urllib2
import numpy
from .. import ChannelConverter
from .. import StreamConverter
from ..TimeseriesFactory import TimeseriesFactory
from ..TimeseriesFactoryException import TimeseriesFactoryException
from ..edge import ObservatoryMetadata

from IAFParser import IAFParser
from IAFWriter import IAFWriter


# pattern for IAF file names
# example bou10jan.bin
IAF_FILE_PATTERN = '%(OBS)s%(yb)s.bin'


def read_url(url):
    """Open and read url contents.

    Parameters
    ----------
    url : str
        A urllib2 compatible url, such as http:// or file://.

    Returns
    -------
    str
        contents returned by url.

    Raises
    ------
    urllib2.URLError
        if any occurs
    """
    response = urllib2.urlopen(url)
    content = None
    try:
        content = response.read()
    except urllib2.URLError, e:
        print e.reason
        raise
    finally:
        response.close()
    return content


class IAFFactory(TimeseriesFactory):
    """TimeseriesFactory for IAF formatted files.

    Parameters
    ----------
    urlTemplate : str
        A string that contains any of the following replacement patterns:
        - '%(OBS)s' uppercase observatory code
        - '%(ym)s' time formatted as YYMMM

    Notes
    -----
    IAF output writes minutes, daily, hourly, and tri-hourly data.
        The Tri-hourly data is the K index.
        The parser reads in all the data, but since geomag-algorithms
        doesn't currently deal with multiple time intervals in the
        same dataset the writer calculates daily and hourly data
        from the minute data.
    --------
    IAFParser
    IAFWriter
    """

    def __init__(self, urlTemplate, observatory=None, channels=None, type=None,
            interval=None, version1_flag=None):
        TimeseriesFactory.__init__(self, observatory, channels, type, interval)
        self.urlTemplate = urlTemplate
        self._version1_flag = version1_flag

    def get_timeseries(self, starttime, endtime, observatory=None,
            channels=None, type='definitive', interval='minute'):
        """Get timeseries data

        Parameters
        ----------
        observatory : str
            observatory code.
        starttime : obspy.core.UTCDateTime
            time of first sample.
        endtime : obspy.core.UTCDateTime
            time of last sample.
        type : {'definitive'}
            data type.
        interval : {'minute'}
            data interval.

        Returns
        -------
        obspy.core.Stream
            timeseries object with requested data.

        Raises
        ------
        TimeseriesFactoryException
            if invalid values are requested, or errors occur while
            retrieving timeseries.
        """
        observatory = observatory or self.observatory
        channels = channels or self.channels
        type = type or self.type
        interval = interval or self.interval
        months = self._get_months(starttime, endtime)
        timeseries = obspy.core.Stream()
        for month in months:
            url = self._get_url(observatory, month, type, interval)
            iafFile = read_url(url)
            timeseries += self.parse_string(iafFile,
                    self._version1_flag, interval)
        # merge channel traces for multiple months
        timeseries.merge()
        # trim to requested start/end time
        timeseries.trim(starttime, endtime)
        return timeseries

    def parse_string(self, iafString, version1_flag, interval='minute'):
        """Parse the contents of a string in the format of an IAF file.

        Parameters
        ----------
        iafString : str
            string containing iaf binary content.
        version1_flag : int
            Flag indicating file is version1
        interval : {'minute', 'second'}
            data interval.

        Returns
        -------
        obspy.core.Stream
            parsed data.
        """
        parser = IAFParser()

        parser.parse(iafString, interval)
        metadata = parser.metadata
        starttime = parser.starttime
        endtime = parser.endtime - 60
        data = parser.data
        stream = obspy.core.Stream()
        for channel in data.keys():
            length = len(data[channel])
            rate = (length - 1) / (endtime - starttime)
            stats = obspy.core.Stats(metadata)
            stats.starttime = starttime
            stats.sampling_rate = rate
            stats.npts = length
            stats.channel = channel
            if channel == 'D':
                data[channel] = ChannelConverter.get_radians_from_minutes(
                    data[channel])
            if channel == 'K':
                stats.sampling_rate = 0.00009259259259
            stream += obspy.core.Trace(data[channel], stats)

        self._clean_G_F_channels(metadata['version'],
                metadata['channels'], stream)
        self._set_metadata(stream)
        return stream

    def put_timeseries(self, timeseries, starttime=None, endtime=None,
            channels=None, type=None, interval=None):
        """Store timeseries data.

        Parameters
        ----------
        timeseries : obspy.core.Stream
            stream containing traces to store.
        starttime : UTCDateTime
            time of first sample in timeseries to store.
            uses first sample if unspecified.
        endtime : UTCDateTime
            time of last sample in timeseries to store.
            uses last sample if unspecified.
        channels : array_like
            list of channels to store, optional.
            uses default if unspecified.
        type : {'definitive', 'provisional', 'quasi-definitive', 'variation'}
            data type, optional.
            uses default if unspecified.
        interval : {'daily', 'hourly', 'minute', 'monthly', 'second'}
            data interval, optional.
            uses default if unspecified.
        """
        if not self.urlTemplate.startswith('file://'):
            raise TimeseriesFactoryException('Only file urls are supported')
        channels = channels or self.channels
        type = type or self.type
        interval = interval or self.interval
        stats = timeseries[0].stats
        observatory = stats.station
        starttime = starttime or stats.starttime
        endtime = endtime or stats.endtime
        months = self._get_months(starttime, endtime)
        for month in months:
            month_filename = self._get_file_from_url(
                    self._get_url(observatory, month, type, interval))
            month_timeseries = self._get_slice(timeseries, month, interval)
            print month_timeseries
            with open(month_filename, 'wb') as fh:
                self.write_file(fh, month_timeseries, channels)

    def write_file(self, fh, timeseries, channels):
        """writes timeseries data to the given file object.

        Parameters
        ----------
        fh: file object
        timeseries : obspy.core.Stream
            stream containing traces to store.
        channels : array_like
            list of channels to store
        """
        self._set_metadata(timeseries)
        IAFWriter().write(fh, timeseries, channels)

    def _clean_G_F_channels(self, version, channels, stream):
        """
        Parameters
        ----------
        version : int
            The version of the iaf file
        channels : str
            The channels from the orientation field of the file.
        stream : obspy.core.StreamConverter
            The timeseries containing the data
        Notes
        -----
        version 1.0: channel 4, is either Fs, or Fv. If its' Fv
            delete Fv.
        version 1.1: channel 4 is either empty or Fs, nothing needs to be done.
        version 2.0: channel 4 is deltaf, but may contain Fs in places where
            Fv can't be calculated, calculate Fs from deltaf and get/remove
            Fs from channel 4 where necessary.
        version 2.1: Same as 2.0, but an orientation of HDZ or XYZ indicates
            no Fs.
        """
        if (version >= 1.0 and version < 1.1) and self._version1_flag:
            for tr in stream.select(channel='F'):
                stream.remove(tr)
        elif version >= 2 and (channels == 'XYZG' or channels == 'HDZG'):
            z = stream.select(channel='Z')[0]
            if channels == 'XYZG':
                x = stream.select(channel='X')[0]
                y = stream.select(channel='Y')[0]
                fv = ChannelConverter.get_computed_f_using_squares(x, y, z)
            elif channels == 'HDZG':
                h = stream.select(channel='H')[0]
                fv = ChannelConverter.get_computed_f_using_squares(h, 0, z)

            deltaf = stream.select(channel='G')[0]
            fs = numpy.subtract(fv, deltaf)

            # copy fs values out of msg and change deltaf to nan
            for value in xrange(0, len(fv)):
                if numpy.isnan(deltaf[value]) and \
                        not numpy.isnan(fv[value]):
                    fs[value] = -deltaf[0][value]
                    deltaf[0][value] = numpy.isnan

            stream += obspy.core.Stream((StreamConverter._get_trace(
                    'F', deltaf.stats, fs), ))

    def _get_file_from_url(self, url):
        """Get a file for writing.

        Ensures parent directory exists.

        Parameters
        ----------
        url : str
            Url path to IMFV283

        Returns
        -------
        str
            path to file without file:// prefix

        Raises
        ------
        TimeseriesFactoryException
            if url does not start with file://
        """
        if not url.startswith('file://'):
            raise TimeseriesFactoryException(
                    'Only file urls are supported for writing')
        filename = url.replace('file://', '')
        parent = os.path.dirname(filename)
        if not os.path.exists(parent):
            os.makedirs(parent)
        return filename

    def _get_months(self, starttime, endtime):
        """Get months between (inclusive) starttime and endtime.

        Parameters
        ----------
        starttime : obspy.core.UTCDateTime
            the start time
        endtime : obspy.core.UTCDateTime
            the end time

        Returns
        -------
        array_like
            list of times, one per months, for all months between and including
            ``starttime`` and ``endtime``.

        Raises
        ------
        TimeseriesFactoryException
            if starttime is after endtime
        """
        if starttime > endtime:
            raise TimeseriesFactoryException(
                    'starttime must be before endtime')
        months = []
        month = starttime
        lastday = (endtime.year, endtime.month)
        while True:
            months.append(month)
            if lastday == (month.year, month.month):
                break
            # move to next day
            month = obspy.core.UTCDateTime(month.year, month.month + 1,
                    month.day)
        return months

    def _get_slice(self, timeseries, month, interval):
        """Get the first and last time for a day

        Parameters
        ----------
        timeseries : obspy.core.Stream
            timeseries to slice
        day : UTCDateTime
            time in day to slice

        Returns
        -------
        obspy.core.Stream
            sliced stream
        """
        month = month.datetime
        start = obspy.core.UTCDateTime(month.year, month.month, 1, 0, 0, 0)
        if interval == 'minute':
            end = start + 86340.0 * 32
            end = obspy.core.UTCDateTime(end.year, end.month, 1, 0, 0, 0)
            end = end - 1
        return timeseries.slice(start, end)

    def _get_url(self, observatory, date,
                type='definitive', interval='minute'):
        """Get the url for a specified IMFV283 file.

        Replaces patterns (described in class docstring) with values based on
        parameter values.

        Parameters
        ----------
        observatory : str
            observatory code.
        date : obspy.core.UTCDateTime
            month to fetch (year, month)
        type : {'definitive'}
            data type.
        interval : {'minute'}
            data interval.
        """
        return self.urlTemplate % {
                'OBS': observatory.upper(),
                'ym': date.strftime('%y%b').lower()}

    def _set_metadata(self, stream):
        """Set metadata
        Parameters
        ----------
        stream: obspy.core.Stream

        Notes
        -----
        copies any values NOT in stats already from DefaultMetadata
        """
        observatory_metadata = ObservatoryMetadata()
        for st in stream:
            stats = st.stats
            tmp = {}
            station = stats['station'].strip()
            metadata = \
                    observatory_metadata.metadata[station]['metadata']
            for key in metadata:
                if key not in stats:
                    stats[key] = metadata[key]
            interval = 'minute'
            if 'interval_specific' in\
                    observatory_metadata.metadata[station]:
                interval_specific = \
                    observatory_metadata.metadata[station]\
                            ['interval_specific']
                for key in interval_specific[interval]:
                    stats[key] = interval_specific[interval][key]
