"""Controller class for geomag algorithms"""

import TimeseriesUtilities as TUtils
import TimeseriesFactoryException
import copy


class Controller(object):
    """Controller for geomag algorithms.

    Parameters
    ----------
    inputFactory: TimeseriesFactory
        the factory that will read in timeseries data
    outputFactory: TimeseriesFactory
        the factory that will output the timeseries data
    algorithm: Algorithm
        the algorithm(s) that will procees the timeseries data

    Notes
    -----
    Has 2 basic modes.
    Run simply sends all the data in a stream to edge. If a startime/endtime is
        provided, it will send the data from the stream that is within that
        time span.
    Update will update any data that has changed between the source, and
        the target during a given timeframe. It will also attempt to
        recursively backup so it can update all missing data.
    """

    def __init__(self, inputFactory, outputFactory, algorithm):
        self._inputFactory = inputFactory
        self._algorithm = algorithm
        self._outputFactory = outputFactory

    def run(self, options):
        """run controller
        Parameters
        ----------
        options: dictionary
            The dictionary of all the command line arguments. Could in theory
            contain other options passed in by the controller.
        """
        input_channels = self._algorithm.get_input_channels()
        algorithm_start, algorithm_end = self._algorithm.get_input_interval(
                options.starttime, options.endtime)

        timeseries = self._inputFactory.get_timeseries(algorithm_start,
            algorithm_end, channels=input_channels)

        processed = self._algorithm.process(timeseries)
        output_channels = self._algorithm.get_output_channels()

        output_channels = self._get_output_channels(output_channels,
                options.outchannels)

        self._outputFactory.put_timeseries(timeseries=processed,
                starttime=options.starttime, endtime=options.endtime,
                channels=output_channels)

    def run_as_update(self, options):
        """Updates data.
        Parameters
        ----------
        options: dictionary
            The dictionary of all the command line arguments. Could in theory
            contain other options passed in by the controller.

        Notes
        -----
        Finds gaps in the target data, and if there's new data in the input
            source, calls run with the start/end time of a given gap to fill
            in.
        It checks the start of the target data, and if it's missing, and
            there's new data available, it backs up the starttime/endtime,
            and recursively calls itself, to check the previous period, to see
            if new data is available there as well. Calls run for each new
            period, oldest to newest.
        """
        input_channels = self._algorithm.get_input_channels()
        output_channels = self._algorithm.get_output_channels()

        output_channels = self._get_output_channels(output_channels,
                options.outchannels)

        timeseries_source = self._inputFactory.get_timeseries(
                options.starttime, options.endtime, channels=input_channels)
        timeseries_target = self._outputFactory.get_timeseries(
                options.starttime, options.endtime, channels=output_channels)

        source_gaps = TUtils.get_timeseries_gaps(
                timeseries_source, input_channels, options.starttime,
                options.endtime)
        target_gaps = TUtils.get_timeseries_gaps(
                timeseries_target, output_channels, options.starttime,
                options.endtime)
        source_gaps = TUtils.get_merged_gaps(source_gaps, input_channels)
        target_gaps = TUtils.get_merged_gaps(target_gaps, output_channels)
        del timeseries_source
        del timeseries_target

        if ((not len(source_gaps) or
                len(source_gaps) and source_gaps[0][0] != options.starttime)
                and
                len(target_gaps) and target_gaps[0][0] == options.starttime):
            new_options = copy.deepcopy(options)
            new_options.starttime = (options.starttime -
                    (options.endtime - options.starttime))
            new_options.endtime = (options.starttime -
                    TUtils.get_seconds_of_interval(options.interval))
            self.run_as_update(new_options)
        for target_gap in target_gaps:
            if not TUtils.gap_is_new_data(source_gaps, target_gap):
                continue
            new_options = copy.deepcopy(options)
            new_options.starttime = target_gap[0]
            new_options.endtime = target_gap[1]
            self.run(new_options)

    def _get_output_channels(self, algorithm_channels, commandline_channels):
        """get output channels

        Parameters
        ----------
        algorithm_channels: array_like
            list of channels required by the algorithm
        commandline_channels: array_like
            list of channels requested by the user
        Notes
        -----
        We want to return the channels requested by the user, but we require
            that they be in the list of channels for the algorithm.
        """
        if commandline_channels is not None:
            for channel in commandline_channels:
                if channel not in algorithm_channels:
                    raise TimeseriesFactoryException(
                        'Output "%s" Channel not in Algorithm'
                            % channel)
            return commandline_channels
        return algorithm_channels
