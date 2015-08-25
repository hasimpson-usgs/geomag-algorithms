"""Parsing methods for the IAF Format."""


import numpy
import struct
import sys
from datetime import datetime
from obspy.core import UTCDateTime
from ..edge import EdgeFactory

# values that represent missing data points in IAF
EIGHTS = numpy.float64('888888')
NINES = numpy.float64('999999')
K_NINES = numpy.float64('999')
IAFSTRING = '<4s4i4s4s1i4s4s2i4s4s1i1i1440i1440i1440i1440i' + \
        '24i24i24i24i1i1i1i1i8i4i'
REC_LENGTH = 23552
DAY_SECONDS = 86400


class IAFParser(object):
    """IAF parser.

    Based on documentation at:
      http://www.intermagnet.org/data-donnee/formats/iaf-eng.php

    Attributes
    ----------
    headers : dict
        parsed IAF headers.
    metadata : dict
        metadata provided
    channels : array
        parsed channel names.
    times : array
        parsed timeseries times.
    data : dict
        keys are channel names (order listed in ``self.channels``).
        values are ``numpy.array`` of timeseries values, array values are
        ``numpy.nan`` when values are missing.
    """

    def __init__(self):
        """Create a new IAF parser."""
        # header fields
        self.headers = {
            'data_type': 'definitive'
        }
        self.metadata = {
            'network': 'NT'
        }
        # array of channel names
        self.channels = []
        # dictionary of data (channel : numpy.array<float64>)
        self.data = {}

    def parse(self, data, interval):
        """Parse a string containing IAF formatted data.

        Parameters
        ----------
        data : str
            IAF binary formatted file contents.
        interval : {'minute', 'hourly', 'daily'}
            interval of data to be returned
        """
        start = 0
        end = REC_LENGTH
        days = int(len(data)/REC_LENGTH)
        for day in range(0, days):
            start = day * REC_LENGTH
            end = (day + 1) * REC_LENGTH
            record = struct.unpack(IAFSTRING, data[start:end])
            self._parse_header(record)
            if day == 0:
                self._parse_channels()
            self._parse_data(interval, record)

        self._post_process()

    def _parse_header(self, record):
        """Parse header line.

        Parameters
        ----------
        record : str
            a single day of data

        Adds value to ``self.headers``
        """
        self.headers['station'] = record[0].strip()
        self.headers['date'] = record[1]
        self.headers['geodetic_latitude'] = str(record[2] / 1000.0)
        self.headers['geodetic_longitude'] = str(record[3] / 1000.0)
        self.headers['elevation'] = str(record[4])
        self.headers['channels'] = record[5]
        self.headers['agency_name'] = record[6]
        self.headers['d_conversion'] = record[7]
        self.headers['data_quality'] = record[8]
        self.headers['instrumentation'] = record[9]
        self.headers['k_9'] = record[10]
        self.headers['sensor_sampling_rate'] = record[11] / 1000.0
        self.headers['sensor_orientation'] = record[12]
        self.headers['publication_date'] = record[13]
        self.headers['version'] = record[14]

    def _parse_channels(self):
        """Parse data header that contains channel names.

        Adds channel names to ``self.channels``.
        """
        self.channels.append(self.headers['channels'][0])
        self.channels.append(self.headers['channels'][1])
        self.channels.append(self.headers['channels'][2])
        if len(self.headers['channels']) > 3:
            self.channels.append(self.headers['channels'][3])
        else:
            self.channels.append('')

        self.channels.append('K')

    def _parse_data(self, interval, record):
        """parse data

        Parameters
        ----------
        interval : {'minute', 'hourly', 'daily'}
            interval of data to be returned
        record : str
            a single day of data

        Adds values to ``self.data``
        """
        if 'station' not in self.metadata:
            self.starttime = UTCDateTime(str(self.headers['date']))
            self.metadata.update(self.headers)
            for channel in self.channels:
                self.data[channel] = []
        self.endtime = UTCDateTime(str(self.headers['date'])) + DAY_SECONDS

        if interval == 'minute':
            for channel, count in zip(self.channels[0:4], xrange(0, 4)):
                start_record = 16 + 1440 * count
                end_record = 16 + 1440 * (count+1)
                self.data[channel].extend(record[start_record:end_record])
        elif interval == 'hourly':
            for channel, count in zip(self.channels[0:4], xrange(0, 4)):
                start_record = 5776 + 24 * (count)
                end_record = 5776 + 24 * (count + 1)
                self.data[channel].extend(record[start_record:end_record])
        elif interval == 'daily':
            for channel, count in zip(self.channels[0:4], xrange(0, 4)):
                start_record = 5872 + count
                end_record = 5872 + (count + 1)
                self.data[channel].extend(record[start_record:end_record])

        # K values
        channel = self.channels[4]
        start_record = 5876
        end_record = 5884
        self.data[channel].extend(record[start_record:end_record])

    def _post_process(self):
        """
        NOTES
        -----
        1) convert all arrays to numpy arrays
        2) convert missing data to nan's
        3) convert all channels but K to correct value by dividing by 10
        """

        for channel in self.data.keys():
            data = numpy.array(self.data[channel], dtype=numpy.float64)
            # filter empty values
            data[data == EIGHTS] = numpy.nan
            data[data == NINES] = numpy.nan
            if channel == 'K':
                data[data == K_NINES] = numpy.nan
            else:
                data = numpy.divide(data, 10.)
            self.data[channel] = data
