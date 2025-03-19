""" Global Ionosphere Map (GIM) Module.
This module implements classes to store and manage (interpolate) GIM data from IONEX files.
"""

import numpy as np
import os
import math

from src.common_log import get_logger, IO_LOG
from src.data_mng import TimeSeries, Container
from src.data_types.date import Epoch
from ...constants import PI, SECONDS_IN_DAY, RAD2DEG
from ...io.config import config_dict
from ...io.rinex_parser.ionex_reader import IONEXReader
from ...utils.interpolation import binary_search


class GIMHeader(Container):
    """
    GIMHeader class, inherits from the :py:class:`Container` class.
    Container class to store the header information of an IONEX file.

    Attributes:
        first_epoch(src.data_types.date.Epoch): first epoch of the GIM
        last_epoch(src.data_types.date.Epoch): last epoch of the GIM
        n_maps(int): number of maps in the GIM file
        mapping_function(str): mapping function used in the GIM file
        base_radius(float): base radius of the Earth in km
        hgt(float): height grid array
        lat(list): latitude grid array
        lon(list): longitude grid array
        exponent(int): exponent of the TEC values
    """

    __slots__ = ["first_epoch", "last_epoch", "n_maps", "mapping_function", "base_radius",
                 "hgt", "lat", "lon", "exponent"]

    def __init__(self):
        super().__init__()
        self.first_epoch = None
        self.last_epoch = None
        self.n_maps = 0
        self.mapping_function = None
        self.base_radius = 0
        self.hgt = None
        self.lat = None
        self.lon = None
        self.exponent = 0

    def __str__(self):
        return f"First epoch: {self.first_epoch}\n" \
               f"Last epoch: {self.last_epoch}\n" \
               f"Number of maps: {self.n_maps}\n" \
               f"Mapping function: {self.mapping_function}\n" \
               f"Base radius: {self.base_radius}\n" \
               f"Height: {self.hgt}\n" \
               f"Latitude: {self.lat}\n" \
               f"Longitude: {self.lon}\n" \
               f"Exponent: {self.exponent}"

    def __repr__(self):
        return self.__str__()


class GlobalIonoMap:
    """
    GlobalIonoMap class to store and manage (interpolate) GIM data from IONEX files.
    """

    def __init__(self):
        self.header = GIMHeader()
        self.lat_array = []
        self.lon_array = []
        self.height = 0
        self.gim_data = TimeSeries()
        self.write_trace = False
        self.trace_file = None

    def __del__(self):
        if self.trace_file is not None:
            self.trace_file.close()

    def __str__(self):
        """ Print the GIM data to a string, for debug purposes """
        metadata = f"{self.header}\n" \
                   f"Latitude array: {self.lat_array}\n" \
                   f"Longitude array: {self.lon_array}\n" \
                   f"Height: {self.height}\n"

        data = ""
        vEpochs = self.gim_data.get_all_epochs()
        for epoch in vEpochs:
            tec_mat = self.gim_data.get_data_for_epoch(epoch)
            with np.printoptions(threshold=np.inf):
                data += f"{epoch}\n{tec_mat}\n"
            data += "--------------------------------\n"
        return metadata + data

    def init(self, ionex_files, trace_dir):
        """
        Initialize this object with the Global Ionosphere Maps (GIM) data from precise IONEX files.

        NOTE: Currently only one IONEX file is supported. The first one will be used.

        Args:
            ionex_files(list): list of IONEX clock files from the user configuration
            trace_dir(str): directory to store the trace file
        """
        log = get_logger(IO_LOG)
        log.info("Launching IONEX GIM Reader...")
        if len(ionex_files) != 1:
            log.warning("Currently only one IONEX file is supported. The first one will be used.")
        IONEXReader(ionex_files[0], self)

        self.write_trace = config_dict.get("solver", "trace_files")
        if self.write_trace:
            inputs_dir = f"{trace_dir}\\interpolations"
            if not os.path.isdir(inputs_dir):
                try:
                    os.makedirs(inputs_dir)
                except:
                    raise IOError(f"Cannot create dir: {inputs_dir}")
            self.trace_file = open(f"{inputs_dir}\\gim_interpolations.txt", "w")

    def initialize_arrays(self):
        """ Initialize internal arrays from the header information. """
        lat_interval = abs(self.header.lat[1] - self.header.lat[0])
        n_lat = int(lat_interval / abs(self.header.lat[2])) + 1
        lon_interval = abs(self.header.lon[1] - self.header.lon[0])
        n_lon = int(lon_interval / abs(self.header.lon[2])) + 1
        self.lat_array = np.linspace(self.header.lat[0], self.header.lat[1], n_lat)
        self.lon_array = np.linspace(self.header.lon[0], self.header.lon[1], n_lon)
        self.height = self.header.hgt[0]

    def set_data(self, epoch, data):
        """
        Method to set a GIM map for the provided epoch.

        Args:
            epoch (Epoch): epoch to set the GIM map
            data (numpy.ndarray): numpy matrix with GIM data for the provided epoch

        Raises:
            AttributeError: this exception is raised if one of the input arguments is not of the valid type.
        """
        if not isinstance(epoch, Epoch):
            raise AttributeError(f'First argument should be a valid Epoch object. Type {type(epoch)} was provided '
                                 f'instead')
        if not isinstance(data, np.ndarray):
            raise AttributeError(f'Second argument should be a valid numpy ndarray. Type {type(data)} '
                                 f'was provided instead')

        self.gim_data[epoch] = data

    def interpolate(self, epoch, lat, lon, sun_fixed=False):
        """
        Method to interpolate the Total Electron Content (TEC) for a given epoch and pierce point coordinates.

        Args:
            epoch (Epoch): epoch to compute the TEC
            lat (float): latitude of the pierce point in [rad]
            lon (float): longitude of the pierce point in [rad]
            sun_fixed(bool): flag to indicate if the longitude should be fixed
                for earth rotation correction (sun-fixed coordinate)

        Returns:
            float: computed VTEC value for the provided epoch and pierce point coordinates
        """
        if self.write_trace:
            self.trace_file.write(f"Computing TEC for epoch {epoch} at pierce point coordinates: "
                                  f"lat {lat * RAD2DEG} [deg], lon {lon * RAD2DEG} [deg].\n")

        # get the two epochs closest to the desired epoch
        vEpochs = self.gim_data.get_n_items(epoch, 1)
        epoch1 = vEpochs[0]
        epoch2 = vEpochs[1]

        tau = (epoch - epoch1).total_seconds() / (epoch2 - epoch1).total_seconds()
        tec_matrix1 = self.gim_data.get_data_for_epoch(epoch1)
        tec_matrix2 = self.gim_data.get_data_for_epoch(epoch2)

        if sun_fixed:
            # update longitude for earth rotation correction (sun-fixed coordinate)
            dt1 = (epoch - epoch1).total_seconds()
            dt2 = (epoch - epoch2).total_seconds()

            lon_epoch1 = lon + 2.0 * PI * dt1 / SECONDS_IN_DAY
            lon_epoch2 = lon + 2.0 * PI * dt2 / SECONDS_IN_DAY

            # make sure that lon_rotated is in the range [-pi,pi] radians
            lon_epoch1 = math.remainder(lon_epoch1, math.tau)
            lon_epoch2 = math.remainder(lon_epoch2, math.tau)
        else:
            lon_epoch1 = math.remainder(lon, math.tau)
            lon_epoch2 = math.remainder(lon, math.tau)

        vtec1 = self._interpolate_for_epoch(lat, lon_epoch1, tec_matrix1)
        vtec2 = self._interpolate_for_epoch(lat, lon_epoch2, tec_matrix2)
        vtec = (1 - tau) * vtec1 + tau * vtec2

        if self.write_trace:
            self.trace_file.write(f"\tComputed VTEC = {vtec}.\n")
        return vtec

    def _interpolate_for_epoch(self, lat_pp: float, lon_pp: float, tec_matrix: np.ndarray) -> float:
        # convert lat_pp and lon_pp to degrees
        lat_pp = lat_pp * RAD2DEG
        lon_pp = lon_pp * RAD2DEG

        # find the closest latitude and longitude grid points and indexes in array
        grid_lat, idx_lat = binary_search(list(self.lat_array), lat_pp, 1, ret_index=True, extrapolation=True)
        grid_lon, idx_lon = binary_search(list(self.lon_array), lon_pp, 1, ret_index=True, extrapolation=True)

        if self.write_trace:
            self.trace_file.write(f"\t[GRID POINT] VTEC for lat {grid_lat[0]} and lon {grid_lon[0]}: "
                                  f"{tec_matrix[idx_lat[0], idx_lon[0]]}\n")
            self.trace_file.write(f"\t[GRID POINT] VTEC for lat {grid_lat[1]} and lon {grid_lon[1]}: "
                                  f"{tec_matrix[idx_lat[1], idx_lon[1]]}\n")
            self.trace_file.write(f"\t[GRID POINT] VTEC for lat {grid_lat[0]} and lon {grid_lon[1]}: "
                                  f"{tec_matrix[idx_lat[0], idx_lon[1]]}\n")
            self.trace_file.write(f"\t[GRID POINT] VTEC for lat {grid_lat[1]} and lon {grid_lon[0]}: "
                                  f"{tec_matrix[idx_lat[1], idx_lon[0]]}\n")

        lat_coeff = (lat_pp - grid_lat[0]) / (grid_lat[1] - grid_lat[0])
        lon_coeff = (lon_pp - grid_lon[0]) / (grid_lon[1] - grid_lon[0])

        if 0 <= lat_coeff <= 1 and 0 <= lon_coeff <= 1:
            vtec = (1 - lat_coeff) * (1 - lon_coeff) * tec_matrix[idx_lat[0], idx_lon[0]] + \
                   lat_coeff * (1 - lon_coeff) * tec_matrix[idx_lat[1], idx_lon[0]] + \
                   (1 - lat_coeff) * lon_coeff * tec_matrix[idx_lat[0], idx_lon[1]] + \
                   lat_coeff * lon_coeff * tec_matrix[idx_lat[1], idx_lon[1]]
            if self.write_trace:
                self.trace_file.write(f"\t[INTERPOLATED] VTEC for lat {lat_pp} and lon {lon_pp}: {vtec}\n")
        else:
            # nearest neighbor extrapolation
            if lat_coeff <= 0.5 and lon_coeff <= 0.5:
                vtec = tec_matrix[idx_lat[0], idx_lon[0]]
            elif lat_coeff <= 0.5 and lon_coeff > 0.5:
                vtec = tec_matrix[idx_lat[0], idx_lon[1]]
            elif lat_coeff > 0.5 and lon_coeff <= 0.5:
                vtec = tec_matrix[idx_lat[1], idx_lon[0]]
            elif lat_coeff > 0.5 and lon_coeff > 0.5:
                vtec = tec_matrix[idx_lat[1], idx_lon[1]]
            else:
                raise ValueError(f"Error in interpolation coefficients: lat_coeff={lat_coeff}, lon_coeff={lon_coeff}")
            if self.write_trace:
                self.trace_file.write(f"\t[EXTRAPOLATED] VTEC for lat {lat_pp} and lon {lon_pp}: {vtec} "
                                      f"(pulled to the nearest neighbor lat_coeff={lat_coeff}, "
                                      f"lon_coeff={lon_coeff})\n")

        return vtec
