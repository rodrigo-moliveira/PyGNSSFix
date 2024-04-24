import numpy as np






class EphemeridePropagator:

    @staticmethod
    def compute_sat_nav_position_dt_rel(nav_message, time_emission, transit) -> \
            tuple[np.array, np.array, float]:
        """
        Computes:
            * the satellite ephemeride
            * the clock relativistic correction
            * the true range rho
        at the requested epoch, and given the closest (valid) navigation data point and transit time
        the corresponding

        Args:
            nav_message (src.data_types.containers.NavigationData.NavigationPointGPS) : navigation data point object
            time_emission (src.data_types.basics.Epoch.Epoch) : Signal emission time (wrt GPS time system)
            transit (float) : Computed transit time in seconds. Used to rotate the computed satellite ephemeride to the
                            ECEF frame at reception time

        Returns:
            tuple [Position, ~numpy.array, float] : returns the computed satellite position at the request epoch, as
                                                    well as the corresponding clock relativistic corrections

        """
        # satellite coordinates in ECEF frame defined at TX time, relativistic correction for satellite clock
        r_sat, v_sat, dt_relative = EphemeridePropagator.compute_nav_sat_eph(nav_message, time_emission)

        # rotation matrix from ECEF TX to ECEF RX (taking into consideration the signal transmission time)
        _R = dcm_e_i(-transit)

        # get satellite position vector at ECEF frame defined at RX time (to be compared with receiver position)
        p_sat = _R @ r_sat

        # get the inertial velocity vector at ECEF frame defined at RX time
        # NOTE: v_sat is the ECEF velocity defined at TX time, expressed in ECEF TX frame.
        # _R @ v_sat rotates this velocity to the ECEF RX frame
        # np.cross(constants.EARTH_ANGULAR_RATE, _R @ r_sat) is the Earth velocity at TX time,
        # defined in the ECEF RX frame, that allows to convert the satellite velocity as required
        # See Eq. (21.29) of Handbook
        v_sat = _R @ v_sat + np.cross(constants.EARTH_ANGULAR_RATE, p_sat)

        return p_sat, v_sat, dt_relative



