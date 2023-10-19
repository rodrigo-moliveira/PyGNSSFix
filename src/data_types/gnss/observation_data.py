from collections import OrderedDict

from src.data_types.gnss.data_type import DataType
from src.data_types.date.date import Epoch
from src.data_types.gnss.satellite import Satellite
from src.data_mng.timeseries import TimeSeries
from src.data_mng.container import Container
from src.errors import NonExistentObservable, EmptyObservationData, DuplicateObservable
from src.data_types.gnss.observation import Observation


class ObservationHeader(Container):
    """
    ObservationHeader class, derived from Container
    stores relevant data from the header section of a rinex gnss_obs file
    """
    __slots__ = ["rinex_version", "satellite_system", "time_system",
                 "receiver_position", "first_epoch", "last_epoch"]

    def __init__(self):
        super().__init__()
        for attr in self.__slots__:
            setattr(self, attr, None)


class EpochData:
    """
    Class EpochData
    Stores all observations for a single epoch. Consists of an OrderedDict:
        * keys -> satellites
        * items -> List of Observation objects with the observations for the given epoch and satellite
    Attributes
        ----------
        _data : OrderedDict
            The DataType (gnss_obs type) of this gnss_obs
    """

    def __init__(self):
        self._data = OrderedDict()

    def set_observable(self, satellite: Satellite, observation: Observation):
        """
        method to set a new gnss_obs

        Args:
            satellite (Satellite)
            observation (Observation)
        Return:
            bool : true if gnss_obs is successfully set. False when this datatype has already been set for
                the given satellite and epoch (overwriting is not legal!)
        """
        if satellite in self._data:
            if observation not in self._data[satellite]:
                self._data[satellite].append(observation)
            else:
                raise DuplicateObservable(
                    f"Trying to set an observable of type {str(observation)} for satellite {str(satellite)},"
                    f"which has already been set. Overwriting not permitted. Ignoring the second appearance")
        else:
            self._data[satellite] = [observation]
        return True

    def get_observables(self, sat: Satellite):
        """
        Args:
                sat (Satellite)
        Return:
                list : list of all observables for the provided satellite (only one epoch)
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
            Observation : gets the gnss_obs for the provided datatype
        Raises:
            NonExistentObservable
        """
        obs_list = self._data[sat]
        for _obs in obs_list:
            if obs == _obs.datatype:
                return _obs
        raise NonExistentObservable(f"gnss_obs {str(obs)} not found for satellite {str(sat)}")

    def has_observable(self, sat, obs):
        obs_list = self._data[sat]
        for _obs in obs_list:
            if obs == _obs.datatype:
                return True
        return False

    def get_satellites(self):
        return list(self._data.keys())

    def get_sats_for_datatypes(self, datatype_list):
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
        if sat in self._data:
            self._data.pop(sat)

    def copy(self):
        # self._data = OrderedDict()
        obj = EpochData()
        for sat, obs_list in self._data.items():
            # construct new gnss_obs and set it
            for obs in obs_list:
                obj.set_observable(sat, obs.copy())
        return obj


class ObservationData:
    """
    Class ObservationData
    Container that stores all batch observations to be processed in the GNSS positioning algorithm.
    Consists of an Ordered Dict:
        * keys -> epochs (Epoch objects)
        * items -> EpochData objects
    Attributes
        ----------
        _data : TimeSeries
        """

    def __init__(self):
        self._data = TimeSeries()
        self._types = {"GPS": [], "GAL": []}
        self._satellites = []
        self.header = ObservationHeader()

    def __str__(self):
        data_str = f"{repr(self.header)}\n{str(self._data)}"

        return data_str

    def set_observable(self, epoch: Epoch, satellite: Satellite, obs_type: DataType, value: float):
        """
        method to set a new gnss_obs (read directly from the rinex gnss_obs file)

        Args:
            epoch (Epoch) : time at reception of signal (time tag from rinex)
            satellite (Satellite)
            obs_type (DataType) : the datatype of the gnss_obs
            value (float) : numeric value of the gnss_obs
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

        # construct gnss_obs
        obs = Observation(obs_type, value)

        self.set_observation(epoch, satellite, obs)

    def set_observation(self, epoch: Epoch, satellite: Satellite, obs: Observation):
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

    def has_type(self, datatype):
        return datatype in self._types[datatype.constellation]

    def has_satellite(self, satellite):
        return satellite in self._satellites

    def remove_observable(self, sat: Satellite, epoch: Epoch, datatype: DataType):
        """
        Remove this observable, for the selected epoch and satellite

         Example: EpochData = [C1, L1, S1, C2, L2, S2]
                remove_observable(datatype = C1)

                -> EpochData = [L1, S1, C2, L2, S2]
        """
        try:
            epoch_data = self._data[epoch]
            epoch_data.remove_observable(sat, datatype)

            # check if there are no satellites (-> remove this epoch_data object)
            if len(epoch_data.get_satellites()) == 0:
                self._data.remove_data(epoch)

        except NonExistentObservable:
            pass

    def remove_for_frequency(self, sat: Satellite, epoch: Epoch, datatype: DataType):
        """
        Remove observations for the selected epoch and satellite which are associated to the frequency
         of the provided datatype.

         Example: EpochData = [C1, L1, S1, C2, L2, S2]
                remove_for_frequency(datatype = C1)

                -> EpochData = [C2, L2, S2]
        """
        try:
            epoch_data = self._data[epoch]
            epoch_data.remove_for_frequency(sat, datatype)
        except NonExistentObservable:
            pass

    # getters
    def get_epoch_data(self, epoch: Epoch):
        """
        Fetch the EpochData object for this epoch

        Args:
            epoch (Epoch)
        Return:
            EpochData  : epoch data object for this epoch
        Raise:
            NonExistentObservable : if the observations are not found
        """
        try:
            return self._data[epoch]
        except KeyError:
            raise NonExistentObservable(f"Non Existent observations for epoch {str(epoch)}")

    def get_observables_at_epoch(self, epoch: Epoch, sat: Satellite):
        """
        Fetch a list of observations for the requested satellite and epoch

        Args:
            sat (Satellite)
            epoch (Epoch)
        Return:
            list : list of observables for the provided sat and epoch
        Raise:
            NonExistentObservable : if the observations are not found
        """

        try:
            return self._data[epoch].get_observables(sat)
        except KeyError:
            raise NonExistentObservable(f"Non Existent gnss_obs for satellite {str(sat)} "
                                        f"and epoch {str(epoch)}")

    def get_observable_at_epoch(self, sat: Satellite, epoch: Epoch, obs: DataType):
        """
        Fetch the requested gnss_obs from the database

        Args:
            sat (Satellite)
            epoch (Epoch)
            obs (DataType)
        Return:
            Observation : the requested gnss_obs
        Raises:
            NonExistentObservable : if the gnss_obs is not found
        """
        try:
            return self._data[epoch].get_observable(sat, obs)
        except KeyError:
            raise NonExistentObservable(f"Non Existent gnss_obs for type {str(obs)}, satellite {str(sat)} "
                                        f"and epoch {str(epoch)}")

    def get_epochs(self):
        return self._data.get_all_epochs()

    def get_satellites(self):
        return self._satellites

    def get_types(self, constellation):
        return self._types[constellation]

    def get_code_types(self, constellation):
        types = list(self.get_types(constellation))
        for t in types:
            if not DataType.is_code(t):
                types.remove(t)
        return types

    def get_satellite_list(self):
        return [str(sat) for sat in self._satellites]

    def get_rate(self):
        epochs = self.get_epochs()
        if len(epochs) > 2:
            return (epochs[1] - epochs[0]).total_seconds()
        raise EmptyObservationData(f"Observation Data is empty")

    def get_first_arc_epoch(self, sat, epoch, rate):
        """Return the first epoch of the arc, given the provided rate
        """
        while True:
            try:
                self.get_observables_at_epoch(epoch, sat)
                epoch = epoch + (-rate)

            except NonExistentObservable:
                return epoch + rate

    def copy(self):
        """ return a copy of this object"""
        obj = ObservationData()

        # shallow/reference copies and deep copy of data
        obj._types = self._types
        obj._satellites = self._satellites
        obj.header = self.header
        obj._data = self._data.copy()

        return obj

    def set_observations_for_constellation(self, constellation, obs_data_in):
        # set observations given an ObservationData object
        epochs = obs_data_in.get_epochs()
        for epoch in epochs:
            epoch_data = obs_data_in.get_epoch_data(epoch)
            sats = epoch_data.get_satellites()
            for sat in sats:
                if sat.sat_system == constellation:
                    observables = epoch_data.get_observables(sat)
                    for obs in observables:
                        self.set_observation(epoch, sat, obs)

    def get_constellations(self):
        consts = []
        for sat in self.get_satellites():
            if sat.sat_system not in consts:
                consts.append(sat.sat_system)
        return consts
