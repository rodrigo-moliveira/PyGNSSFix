from PositioningSolver.src.gnss.observation_models.atmosphere_obs import ionosphereCorrection, troposphericCorrection
from PositioningSolver.src.gnss.observation_models.clock_obs import SVBroadcastCorrection
from PositioningSolver.src.data_types.basics.DataType import DataType, DataTypeFactory
from PositioningSolver.src.gnss.data_types.Observation import Observation
from PositioningSolver.src.math_utils.Constants import Constant

# Disclaimer: The most correct way to implement the 'ObservationReconstruction' class would be the following. First,
# define new classes, one for each contribution in the GNSS observation equation (true range, clock, atmosphere, etc.),
# initialize this classes with the appropriate data and define a compute() method for each, which would compute the
# corresponding contribution to the observation equation.
# Then, in this ObservationReconstruction, store the list of models in the constructor. Then we would simply need to
# iterate over all available models and call the compute() method. By far, this is the cleanest way.
# Since, I only have one model (the GPS recommended one) for each contribution, I can get away with boolean variables,
# which tell us whether or not to include the following contribution.

C1 = DataTypeFactory("C1")
f1 = C1.freq


class ObservationReconstruction:
    def __init__(self, system_geometry, datatype, tropo=False, iono=False, relativistic_correction=False,
                 true_range=None, satellite_clock=False):
        self._model = {"tropo": tropo,  # a priori troposphere model
                       "iono": iono,  # a priori ionosphere model
                       "relativistic_correction": relativistic_correction,
                       "true_range": true_range,
                       "satellite_clock": satellite_clock
                       }
        self._datatype = datatype
        self._system_geometry = system_geometry

    def compute(self, nav_message, nav_header, sat, epoch):
        obs = 0

        # true range
        if self._model["true_range"]:
            obs += self._system_geometry.get("true_range", sat)
        # satellite clock
        if self._model["satellite_clock"]:

            # get TGD and fix it for non L1 users
            TGD = nav_message.TGD
            if DataType.is_iono_free_smooth_code(self._datatype) or DataType.is_iono_free_code(self._datatype):
                TGD = 0  # TGD is 0 for Iono Free observables
            elif self._datatype.freq != f1:  # correct TGD for L2 users
                TGD = (f1.freq_value / self._datatype.freq.freq_value) ** 2 * TGD

            dt_sat, _ = SVBroadcastCorrection(nav_message.af0,
                                              nav_message.af1,
                                              nav_message.af2,
                                              nav_message.toc,
                                              self._system_geometry.get("time_emission", sat))

            if self._model["relativistic_correction"]:
                dt_sat += self._system_geometry.get("dt_rel_correction", sat)

            # correct for TGD (already corrected for the appropriate frequency, and is 0 for IF observables)
            obs -= (dt_sat - TGD) * Constant.SPEED_OF_LIGHT  # convert dt_sat from seconds to meters using c

        # ionosphere
        if self._model["iono"]:
            receiver_position = self._system_geometry.get("receiver_position", sat)
            receiver_position.form = "geodetic"
            az = self._system_geometry.get("az", sat)
            el = self._system_geometry.get("el", sat)
            time_reception = self._system_geometry.get("time_reception", sat)

            obs += ionosphereCorrection(receiver_position[0], receiver_position[1], el, az,
                                        nav_header.iono_corrections["GPSA"],
                                        nav_header.iono_corrections["GPSB"],
                                        time_reception, frequency=self._datatype.freq)

        # troposphere
        if self._model["tropo"]:
            receiver_position = self._system_geometry.get("receiver_position", sat)
            receiver_position.form = "geodetic"
            el = self._system_geometry.get("el", sat)

            obs += troposphericCorrection(receiver_position[2], receiver_position[0], epoch.to_DOY(), el)

        return Observation(self._datatype, obs)
