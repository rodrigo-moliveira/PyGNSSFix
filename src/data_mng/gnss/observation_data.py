from collections import OrderedDict
from datetime import timedelta

from src.data_types.gnss import DataType, Satellite, Observation
from src.data_types.date import Epoch
from src.data_mng import TimeSeries
from src.errors import TimeSeriesError

__all__ = ["ObservationData"]


class EpochData:
    """
    Stores all observations for a single epoch. Consists of an OrderedDict:
        * keys -> satellites
        * Values -> List of Observation objects with the observations for the given epoch and satellite
    """

    def __init__(self):
        self._data = OrderedDict()

    def set_observable(self, satellite: Satellite, observation: Observation):
        """
        set a new observation for the provided satellite

        Args:
            satellite (Satellite)
            observation (Observation)

        Raises:
            TimeSeriesError: raises an exception if the observation has already been set (overwriting is not legal!)
        """
        if satellite in self._data:
            if observation not in self._data[satellite]:
                self._data[satellite].append(observation)
            else:
                raise TimeSeriesError(
                    f"Trying to set an observable of type {str(observation)} for satellite {str(satellite)},"
                    f"which has already been set. Overwriting not permitted. Ignoring the second appearance")
        else:
            self._data[satellite] = [observation]

    def get_observables(self, sat: Satellite):
        """
        Return a list with all observations for the provided satellite
        Args:
            sat (Satellite)
        Return:
            list: list with required observations
        """
        return self._data[sat]

    def get_code_observables(self, sat: Satellite):
        """
        Args:
            sat (Satellite)
        Return:
            list : list of all code (pseudorange) observables for the provided satellite (only one epoch)
        """
        observables = []
        for obs in self._data[sat]:
            if DataType.is_code(obs.datatype):
                observables.append(obs)
        return observables

    def get_observable(self, sat, obs):
        """
        Args:
            sat (Satellite)
            obs (DataType)
        Return:
            Observation : returns the required observation, or raises a TimeSeriesError
        """
        obs_list = self._data[sat]
        for _obs in obs_list:
            if obs == _obs.datatype:
                return _obs
        raise TimeSeriesError(f"Observation {str(obs)} not found for satellite {str(sat)}")

    def has_observable(self, sat: Satellite, obs: Observation):
        """Returns True if the provided Observation is available

        Return:
            bool: True if the observation is available
        """
        obs_list = self._data[sat]
        for _obs in obs_list:
            if obs == _obs.datatype:
                return True
        return False

    def get_satellites(self):
        """
        Return:
            list: returns a list with all the available satellites
        """
        return list(self._data.keys())

    def get_satellites_for_constellation(self, constellation):
        """
        Return:
            list: returns a list with all the available satellites for the provided constellation
        """
        return [sat for sat in list(self._data.keys()) if constellation == sat.sat_system]

    def get_sats_for_datatypes(self, datatype_list):
        """
        Return:
            list: returns a list with all the available satellites for the provided Observations
        """
        sat_list = []

        for sat in self.get_satellites():
            has_type = True  # we assume that this satellite has this type

            # iterate over all requested types
            for datatype in datatype_list:
                if datatype is not None and not self.has_observable(sat, datatype):
                    has_type = False  # we found out that actually this satellites does NOT have this type

            # if this satellite has data for all requested datatypes, append to the list
            if has_type:
                sat_list.append(sat)

        return sat_list

    def __str__(self):
        data_str = ""
        for sat, obs in self._data.items():
            data_str += "\t" + str(sat) + " -> " + str(obs) + "\n"

        if data_str == "":
            data_str = "\t-->Empty Epoch Data"

        return data_str

    def remove_observable(self, sat: Satellite, datatype: DataType):
        """Removes all the observations associated with the provided Satellite and DataType"""
        obs_list = self._data[sat]
        new_obs_list = []

        for obs in obs_list:
            if obs.datatype != datatype:
                new_obs_list.append(obs)
        if new_obs_list:
            self._data[sat] = new_obs_list
        else:
            self._data.pop(sat)

    def remove_for_frequency(self, sat: Satellite, datatype: DataType):
        """Removes all the observations associated with the provided Satellite and frequency"""
        obs_list = self._data[sat]
        new_obs_list = []

        for obs in obs_list:
            if obs.datatype.freq != datatype.freq:
                new_obs_list.append(obs)
        if new_obs_list:
            self._data[sat] = new_obs_list
        else:
            self._data.pop(sat)

    def remove_satellite(self, sat):
        """Removes all the observations associated with the provided Satellite"""
        if sat in self._data:
            self._data.pop(sat)

    def copy(self):
        # self._data = OrderedDict()
        obj = EpochData()
        for sat, obs_list in self._data.items():
            for obs in obs_list:
                obj.set_observable(sat, obs.copy())
        return obj


class ObservationData:
    """
    Observation Data Frame

    Consists of an Ordered Dict:
        * keys -> epochs (Epoch objects)
        * values -> EpochData objects
    """

    def __init__(self):
        self._data = TimeSeries()
        self._types = {"GPS": [], "GAL": []}  # for each key, store all available DataTypes
        self._satellites = []  # list with all available satellites

    def __str__(self):
        data_str = f"{str(self._data)}"

        return data_str

    # Setter Methods
    def set_observable(self, epoch: Epoch, satellite: Satellite, obs_type: DataType, value: float):
        """
        method to set a new Observation

        Args:
            epoch (Epoch) : time at reception of signal (time tag from rinex_parser)
            satellite (Satellite)
            obs_type (DataType) : the datatype of the observation
            value (float) : numeric value of the observation
        """
        if not isinstance(epoch, Epoch):
            raise TypeError(f'First argument should be a valid Epoch object. Type {type(epoch)} was provided instead')
        if not isinstance(satellite, Satellite):
            raise TypeError(f'Second argument should be a valid Satellite object. Type {type(satellite)} '
                            f'was provided instead')
        if not isinstance(obs_type, DataType):
            raise TypeError(f'Third argument should be a valid DataType object. Type {type(obs_type)} '
                            f'was provided instead')
        if not isinstance(value, float) and not isinstance(value, int):
            raise TypeError(f'Forth argument should be a valid number (float or integer). Type {type(value)} '
                            f'was provided instead')

        if isinstance(value, int):
            value = float(value)

        # construct Observation object
        obs = Observation(obs_type, value)

        self.set_observation(epoch, satellite, obs)

    def set_observation(self, epoch: Epoch, satellite: Satellite, obs: Observation):
        """
        method to set a new Observation

        Args:
            epoch (Epoch) : time at reception of signal (time tag from rinex_parser)
            satellite (Satellite)
            obs (Observation) : the Observation object to be set
        """
        if not self._data.has_epoch(epoch):
            # create EpochData object and fill it
            epoch_data = EpochData()
            self._data.set_data(epoch, epoch_data)

        # fetch EpochData object
        epoch_data = self._data[epoch]
        epoch_data.set_observable(satellite, obs)

        # append this obs type to the list
        if obs.datatype not in self._types[satellite.sat_system]:
            self._types[satellite.sat_system].append(obs.datatype)

        if satellite not in self._satellites:
            self._satellites.append(satellite)

    def set_observations_for_constellation(self, constellation, other):
        """
        method to copy all observations from another :py:class:`ObservationData` object
        for the provided constellation

        Args:
            constellation (str) : time at reception of signal (time tag from rinex_parser)
            other (ObservationData): the input ObservationData object
        """
        # set observations given an ObservationData object
        epochs = other.get_epochs()
        for epoch in epochs:
            epoch_data = other.get_epoch_data(epoch)
            sats = epoch_data.get_satellites()
            for sat in sats:
                if sat.sat_system == constellation:
                    observables = epoch_data.get_observables(sat)
                    for obs in observables:
                        self.set_observation(epoch, sat, obs)

    # Getter Methods
    def get_epoch_data(self, epoch: Epoch):
        """
        Fetch the EpochData object for this epoch

        Args:
            epoch (Epoch)
        Return:
            EpochData  : epoch data object for this epoch
        Raises:
            TimeSeriesError : if the observations are not found
        """
        try:
            return self._data[epoch]
        except KeyError:
            raise TimeSeriesError(f"Non Existent observations for epoch {str(epoch)}")

    def get_observables_at_epoch(self, epoch: Epoch, sat: Satellite, datatypes: list = None):
        """
        Fetch a list of observations for the requested satellite and epoch

        Args:
            sat (Satellite)
            epoch (Epoch)
            datatypes (list)
        Return:
            list : list of observables for the provided sat and epoch. If the datatype list is provided, then the return
            list is filtered accordingly
        Raises:
            TimeSeriesError : if the observations are not found
        """

        try:
            out_list = self._data[epoch].get_observables(sat)
            if datatypes is None:
                return out_list
            return [obs for obs in out_list if obs in datatypes]
        except KeyError:
            raise TimeSeriesError(f"Non Existent Observation for satellite {str(sat)} "
                                  f"and epoch {str(epoch)}")

    def get_observable_at_epoch(self, sat: Satellite, epoch: Epoch, obs: DataType):
        """
        Fetch the requested Observation from the database

        Args:
            sat (Satellite)
            epoch (Epoch)
            obs (DataType)
        Return:
            Observation : the requested Observation
        Raises:
            TimeSeriesError : if the Observation is not found
        """
        try:
            return self._data[epoch].get_observable(sat, obs)
        except KeyError:
            raise TimeSeriesError(f"Non Existent Observation for type {str(obs)}, satellite {str(sat)} "
                                  f"and epoch {str(epoch)}")

    def get_epochs(self):
        return self._data.get_all_epochs()

    def get_types(self, constellation):
        return self._types[constellation]

    def get_code_types(self, constellation):
        types = list(self.get_types(constellation))
        code_types = [x for x in types if DataType.is_code(x)]
        return code_types

    def get_doppler_types(self, constellation):
        types = list(self.get_types(constellation))
        code_types = [x for x in types if DataType.is_doppler(x)]
        return code_types

    def get_rate(self):
        epochs = self.get_epochs()
        if len(epochs) > 2:
            return (epochs[1] - epochs[0]).total_seconds()
        raise TimeSeriesError(f"Observation Data is empty")

    def get_first_arc_epoch(self, sat, epoch, rate):
        """Return the first epoch of the arc, given the provided rate
        """
        while True:
            try:
                self.get_observables_at_epoch(epoch, sat)
                epoch = epoch + timedelta(seconds=-rate)

            except TimeSeriesError:
                return epoch + timedelta(seconds=rate)

    # Remove Methods
    def remove_observable(self, sat: Satellite, epoch: Epoch, datatype: DataType):
        """
        Remove the observable, for the selected epoch, satellite and datatype

        Example:
            IN: EpochData = [C1, L1, S1, C2, L2, S2]
            remove_observable(datatype = C1)
            OUT: EpochData = [L1, S1, C2, L2, S2]
        """
        try:
            epoch_data = self._data[epoch]
            epoch_data.remove_observable(sat, datatype)

            # check if there are no satellites (-> remove this epoch_data object)
            if len(epoch_data.get_satellites()) == 0:
                self._data.remove_data(epoch)

        except TimeSeriesError:
            pass

    def remove_for_frequency(self, sat: Satellite, epoch: Epoch, datatype: DataType):
        """
        Remove observations for the selected epoch and satellite which are associated to the frequency
        of the provided datatype.
        """
        try:
            epoch_data = self._data[epoch]
            epoch_data.remove_for_frequency(sat, datatype)
        except TimeSeriesError:
            pass

    # Utilities
    def has_type(self, datatype):
        return datatype in self._types[datatype.constellation]

    def has_satellite(self, satellite):
        return satellite in self._satellites

    def copy(self):
        """ return a copy of this object"""
        obj = ObservationData()

        # shallow/reference copies and deep copy of data
        obj._types = self._types
        obj._satellites = self._satellites
        obj._data = self._data.copy()

        return obj

    def to_csv_file(self):
        """
        Export this ObservationData to a csv file format.
        This function returns the header and body of the file in the following format:
            Week_Number(SCALE),Time_of_Week[s],Satellite,Observation,Measurement,Noise_Std

        Returns:
            str: output string with observation data for this object in csv string format
        """
        if self._data.is_empty():
            return "empty observation data"
        vEpochs = self.get_epochs()
        epoch_scale = vEpochs[0].scale

        # header
        out_string = f"Week_Number({epoch_scale}),Time_of_Week[s],Satellite,Observation,Measurement,Noise_Std\n"

        # body
        for epoch in self.get_epochs():
            week, sow = epoch.gnss_time
            epoch_data = self.get_epoch_data(epoch)
            vSats = epoch_data.get_satellites()
            for sat in vSats:
                vObservables = epoch_data.get_observables(sat)
                for obs in vObservables:
                    out_string = out_string + f"{week},{sow},{sat},{obs.datatype},{obs.value},{obs.std}\n"
        return out_string
