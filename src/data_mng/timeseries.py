from collections import OrderedDict

from src.errors import TimeSeriesError

__all__ = ["TimeSeries"]


class TimeSeries(OrderedDict):
    """
    TimeSeries class, inherits from OrderedDict

    Stores a time series. This class inherits from OrderedDict (dict ordered by insertion order), but
    a feature to sort it by keys is introduced. This is suitable for time series, where keys represent time instants
    (epochs).

    The keys can either be :class:`src.data_types.date.date.Epoch` or floats
    """

    def __init__(self, overwriting=False):
        """
        Construct a Timeseries

        Args:
            overwriting (bool): if True, then data is overwritten when using function :func:`set_data`.
            that has already been set
        """
        super().__init__()
        self._epochs = []
        self._sorted = False
        self.overwriting = overwriting

    # methods to set data
    def set_data(self, epoch, epoch_data):
        """
        Sets the provided data in the TimeSeries frame
        """
        # if not self.overwriting and epoch not in self._epochs:
        if not self.overwriting and epoch in self._epochs:
            # epoch is already defined, and we are in non-overwriting mode -> skip insertion
            raise TimeSeriesError(f"Epoch {epoch} has already been defined. This timeseries is non-overwriting.")

        # add data to internal dict
        self._set_data(epoch, epoch_data)

    def _set_data(self, epoch, epoch_data):
        # add epoch_data for this epoch
        self._epochs.append(epoch)
        self._sorted = False
        super().__setitem__(epoch, epoch_data)

    def __setitem__(self, key, value):
        self.set_data(key, value)

    # method to remove data
    def remove_data(self, epoch):

        if epoch not in self.epochs:
            raise KeyError(f"Key {epoch} not in TimeSeries")

        self.epochs.remove(epoch)
        return self.pop(epoch)

    # getters
    def get_all_epochs(self):
        self._sort()
        return self.epochs

    def get_data_for_epoch(self, epoch):
        self._sort()
        return self[epoch]

    def get_closest_epoch(self, epoch):
        """Returns the closest epoch in the TimeSeries with respect to the provided epoch (argument)"""
        self._sort()

        if len(self.epochs) == 0:
            raise TimeSeriesError(f"TimeSeries is empty.")

        if epoch < self.epochs[0]:
            raise TimeSeriesError(f"Epoch {str(epoch)} is no inside TimeSeries interval "
                                  f"{[repr(i) for i in self.epochs]}")

        # at this point, we can safely assume that vEpochs[0] <= epoch
        prev_epoch = self.epochs[0]
        for this_epoch in self.epochs:
            if this_epoch > epoch:
                return prev_epoch
            prev_epoch = this_epoch

        return prev_epoch

    @staticmethod
    def get_common_epochs(series1, series2):
        """
        Method to align the two time series to a common time interval
        returns a list with the common epochs for both series

        Args:
        series1 (TimeSeries): time series 1
        series2 (TimeSeries): time series 2

        Return:
            list : A list with the common epochs
        """
        # series1.sort()
        # series2.sort()

        # print(series1.epochs[0], series1.epochs[-1])
        # print(series2.epochs[0], series2.epochs[-1])

        epochs = []
        for key in series1.epochs:
            if key in series2.epochs:
                epochs.append(key)

        return epochs

    # utility methods
    def _sort(self):
        # this sort is not very efficient...
        if self.sorted is False:
            epochs = sorted(self.epochs)
            self.epochs.clear()

            new_dct = OrderedDict((key, self[key]) for key in epochs)
            self.clear()
            self.update(new_dct)

            self._sorted = True

    def __repr__(self):
        self._sort()

        my_str = ""
        for epc, obs in self.items():
            my_str += repr(epc) + " -> " + str(obs) + "\n"
        return my_str

    def copy(self):
        """Returns a deep copy of this TimeSeries object"""
        self._sort()
        # clone this TimeSeries
        _copy = TimeSeries()

        # clone data to new object
        for epoch, value in self.items():

            # try to create a clone (deep copy) of the value object (if possible)
            if hasattr(value, "copy"):
                _copy.set_data(epoch, value.copy())
            else:
                _copy.set_data(epoch, value)

        return _copy

    def has_epoch(self, epoch):
        return epoch in self.epochs

    # def get_sub(self, n):
    #    tmOut = []
    #
    #    for i in range(n):
    #        tmSeries = TimeSeries()
    #        for epoch, obs in self._data.items():
    #            tmSeries.setEpochData(epoch, obs[i])
    #        tmOut.append(tmSeries)
    #
    #    return tmOut

    # properties
    @property
    def epochs(self):
        return self._epochs

    @property
    def sorted(self):
        return self._sorted

    def export2time_data(self):
        self._sort()
        return list(self.keys()), list(self.values())

    def is_empty(self):
        return len(self.get_all_epochs()) == 0

    def get_n_items(self, epoch, order):
        """
        Applies binary search to return the `order` number of epochs before and after the provided `epoch`
        This method is helpful to get the surrounding knots to perform interpolations

        Parameters:
            epoch (:class:`src.data_types.date.date.Epoch` or float): The epoch for which to find the surrounding knots.
            order(int): order of the search, i.e., number of epochs before and after to be returned

        Returns:
            tuple: A tuple containing the knot epochs
                   If `epoch` is outside the range of data, a `TimeSeriesError` exception is raised.
        """
        self._sort()

        keys = list(self.keys())

        # Binary search to find the index where `epoch` would be inserted to maintain sorted order
        low, high = 0, len(keys) - 1
        while low <= high:
            mid = (low + high) // 2
            if keys[mid] < epoch:
                low = mid + 1
            else:
                high = mid - 1

        # `low` is now the index where `epoch` would be inserted
        insert_index = low

        # Find the start and end indexes for the sublist
        start_index = max(0, insert_index - order)
        end_index = min(len(keys), insert_index + order)

        # Check if there are enough items before and after `epoch` to return `order` items each way
        if end_index - start_index < 2 * order:
            raise TimeSeriesError(
                f"Not enough elements in the dataset around the specified epoch {epoch} to return the desired number of"
                f" epochs before and after (selected order is {order}).")

        return keys[start_index:insert_index] + keys[insert_index:end_index]
