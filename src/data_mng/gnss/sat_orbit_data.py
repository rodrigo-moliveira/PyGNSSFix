from collections import OrderedDict
import numpy
import numpy as np

from src import constants
from src.data_types.date import Epoch
from src.data_types.gnss import Satellite
from src.data_mng import TimeSeries
from src.models.frames import dcm_e_i
from src.models.gnss_models import compute_nav_sat_eph
from src.io.rinex_parser import SP3OrbitReader
from src.common_log import IO_LOG, get_logger
from src.utils.interpolation import lagrange_interpolation, lagrange_interpolation_derivative

__all__ = ["SatelliteOrbits"]


class SatelliteOrbits:
    """
    Satellite Orbits DataFrame class
    this class stores and manages satellite orbits, either provided from:
        * navigation message (RINEX NAV)
        * precise products (SP3)
    """

    def __init__(self):
        """Constructor of the `SatelliteOrbits` object"""
        self._data = OrderedDict()  # contains the SP3 data
        self.nav_data = None  # contains the broadcast navigation data
        self.use_precise_products = False
        self.interp_order = 0

    def init(self, nav_data, sp3_files, use_precise_products, interp_order=9, first_epoch=None, last_epoch=None):
        """
        Initialize this SatelliteOrbits instance with navigation and precise data.
        The argument `use_precise_products` allows to select which satellite orbits to output.
        If it is set to True, then the precise SP3 data is used, if set to False the
        navigation orbits are used instead

        Args:
            nav_data(src.data_mng.gnss.navigation_data.NavigationData): the `NavigationData` with navigation data
                (must be already initialized)
            sp3_files(list): list of SP3 files from the user configuration
            use_precise_products(bool): True to use precise orbits and False to use broadcast navigation orbits
            interp_order(int): interpolation order for satellite clocks (default order is 9)
            first_epoch(str or None): first observation epoch
            last_epoch(str or None): last observation epoch
        """
        self.nav_data = nav_data
        self.use_precise_products = use_precise_products
        self.interp_order = interp_order

        if use_precise_products:
            log = get_logger(IO_LOG)
            log.info("Launching SP3OrbitReader")
            for file in sp3_files:
                SP3OrbitReader(file, self, first_epoch, last_epoch)

    def __str__(self):
        """Print the SP3 data to a string, for debug purposes"""
        myStr = f"SatelliteOrbits:\n"
        for sat, data in self._data.items():
            myStr += str(sat) + f"\n{str(data)}\n"

        return myStr

    def set_data(self, epoch: Epoch, satellite: Satellite, position: numpy.ndarray):
        """
        method to set a satellite orbit data point for a given epoch and satellite
        Args:
            epoch (Epoch)
            satellite (Satellite)
            position (numpy.ndarray) : satellite orbit from the SP3 Clock file (ECEF cartesian components)
        """

        if not isinstance(epoch, Epoch):
            raise AttributeError(f'First argument should be a valid Epoch object. Type {type(epoch)} was provided '
                                 f'instead')
        if not isinstance(satellite, Satellite):
            raise AttributeError(f'Second argument should be a valid Satellite object. Type {type(satellite)} '
                                 f'was provided instead')
        if not isinstance(position, numpy.ndarray):
            raise AttributeError(f'Third argument should be a numpy.ndarray object. Type {type(position)} '
                                 f'was provided instead')

        if satellite in self._data:
            self._data[satellite].set_data(epoch, position)
        else:
            timeseries = TimeSeries()
            timeseries.set_data(epoch, position)
            self._data[satellite] = timeseries

    # Getters
    def get_data(self):
        return self._data

    def get_sat_data(self, sat):
        try:
            return self._data[sat]
        except KeyError:
            raise KeyError(f"The provided SP3 files do not have data for satellite {str(sat)}")

    def get_satellites(self):
        """Returns a list with all available satellites (from SP3 file)"""
        return list(self._data.keys())

    def get_epochs(self):
        """Returns a list with all available epochs (from SP3 file)"""
        for key, val in self._data.items():
            return val.get_all_epochs()

    def get_orbit(self, sat, epoch):
        """
        Compute the satellite orbits (position, velocity and clock relativistic correction). These quantities are
        computed in the ECEF frame.

        If `use_precise_products` is True then these quantities are computed from the SP3 files. Otherwise,
        they are computed from the provided navigation data.

        Args:
            sat(src.data_types.gnss.Satellite)
            epoch(src.data_types.date.Epoch)
        Returns:
            tuple[numpy.ndarray,numpy.ndarray,float]: tuple with satellite position, velocity and clock relativistic
                correction for the provided satellite and epoch
        """
        if self.use_precise_products:
            return self.get_orbit_precise(sat, epoch)
        else:
            return self.get_orbit_broadcast(sat, epoch)

    def get_orbit_precise(self, sat, epoch):
        """
        Compute the precise satellite orbits from the provided SP3 Orbit data.
        Applies lagrange interpolation when needed (when the provided epoch is in-between data points).

        Args:
            sat(src.data_types.gnss.Satellite)
            epoch(src.data_types.date.Epoch)
        Returns:
            tuple[numpy.ndarray,numpy.ndarray,float]: tuple with satellite position, velocity and clock relativistic
                correction for the provided satellite and epoch

        Raises an exception if the provided satellite does not have valid data
        """
        sat_data = self.get_sat_data(sat)

        n_points = (self.interp_order + 1) // 2
        knots = sat_data.get_n_items(epoch, n_points)

        # check if the epoch is already contained in the TimeSeries or if an interpolation is needed
        if sat_data.has_epoch(epoch):
            pos = sat_data.get_data_for_epoch(epoch).copy()
        else:
            # an interpolation is needed
            poly_pos = lagrange_interpolation([(x-knots[0]).total_seconds() for x in knots],
                                              [sat_data.get_data_for_epoch(x).copy() for x in knots],
                                              self.interp_order)
            pos = poly_pos((epoch-knots[0]).total_seconds())

        # velocity interpolation
        poly_vel = lagrange_interpolation_derivative(
            [(x-knots[0]).total_seconds() for x in knots],
            [sat_data.get_data_for_epoch(x).copy() for x in knots], self.interp_order)
        vel = poly_vel((epoch-knots[0]).total_seconds())

        # computing clock relativistic correction
        rel_correction = -2 * np.dot(pos, vel) / (constants.SPEED_OF_LIGHT**2)

        return pos, vel, rel_correction

    def get_orbit_broadcast(self, sat, epoch):
        """
        Compute the broadcast satellite orbits from the provided RINEX Navigation data.

        Args:
            sat(src.data_types.gnss.Satellite)
            epoch(src.data_types.date.Epoch)
        Returns:
            tuple[numpy.ndarray,numpy.ndarray,float]: tuple with satellite position, velocity and clock relativistic
                correction for the provided satellite and epoch

        Raises an exception if the provided satellite does not have valid data
        """
        nav_message = self.nav_data.get_closest_message(sat, epoch)
        position, velocity, rel_correction = compute_nav_sat_eph(nav_message, epoch)
        # rel_correction2 = -2 * np.dot(position, velocity) / (constants.SPEED_OF_LIGHT ** 2)
        return position, velocity, rel_correction

    def compute_orbit_at_rx_time(self, sat, time_emission, transit):
        """
        This function computes:
            * the satellite position and velocity
            * the clock relativistic correction
        at the requested epoch ('time_emission': time of signal emission). It is noted that the position and velocity
        vectors are rotated to the correct ECEF frame of the reception time, according to the provided transit time.
        That is, the satellite position dnd velocity are computed for the `time_emission` epoch (with respect to the
        ECEF frame defined at this epoch), but then these vectors are rotated to the ECEF frame valid for the reception
        epoch (according to reception_epoch = time_emission + transit).

        This function is called for the reconstruction of the receiver-satellite true range and range rate, since both
        the receiver and satellite vectors must refer to a common ECEF frame

        Args:
            sat (src.data_types.gnss.Satellite)
            time_emission (src.data_types.date.Epoch) : Signal emission time (wrt GPS time system)
            transit (float) : Computed transit time in seconds. Used to rotate the computed satellite ephemeride to the
                            ECEF frame at reception time

        Returns:
            tuple [numpy.ndarray,numpy.ndarray,float] : returns the computed satellite position and velocity vectors
                at the transmission epoch with respect to the ECEF frame of the reception epoch. The satellite
                relativistic clock is also returned.
        """
        # satellite coordinates in ECEF frame defined at TX time, relativistic correction for satellite clock
        r_sat, v_sat, rel_correction = self.get_orbit(sat, time_emission)

        # rotation matrix from ECEF TX to ECEF RX (taking into consideration the signal transmission time)
        _R = dcm_e_i(-transit)

        # get satellite position vector at ECEF frame defined at RX time (to be compared with receiver position)
        p_sat = _R @ r_sat

        # get the inertial velocity vector at ECEF frame defined at RX time
        # NOTE: v_sat is the ECEF velocity defined at TX time, expressed in ECEF TX frame.
        # _R @ v_sat rotates this velocity to the ECEF RX frame
        # np.cross(constants.EARTH_ANGULAR_RATE, _R @ r_sat) is the Earth velocity at TX time,
        # defined in the ECEF RX frame, that allows to convert the satellite velocity as required
        # Reference: Eq. (21.29) of Springer Handbook of Global Navigation Satellite Systems, Springer Cham, 2017
        v_sat = _R @ v_sat + np.cross(constants.EARTH_ANGULAR_RATE, p_sat)

        return p_sat, v_sat, rel_correction
