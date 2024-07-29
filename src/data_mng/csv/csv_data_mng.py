from src.io.states import OUTPUT_FILENAME_MAP
from src.data_mng.container import Container

from .csv_data import CSVData


class GnssRunStorageManager(Container):
    __rx_based__ = ["time", "position", "velocity", "clock_bias",
                    "dop_ecef", "dop_local",
                    "clock_bias_rate", "isb", "tropo_wet"]
    __sat_based__ = ["prefit_residuals", "postfit_residuals",
                     "iono", "satellite_azel", "vel_prefit_residuals", "vel_postfit_residuals"]

    __slots__ = __rx_based__ + __sat_based__ + ["_available"]

    def __init__(self):
        super().__init__()

        # time
        self.time = CSVData(name="time",
                            description="Time in GPS Time Scale",
                            units=None,
                            legend=None,
                            title="GPS Time Scale",
                            time_cols=(0, 1),
                            data_cols=(2,))

        # position in ECEF
        self.position = CSVData(name="position",
                                description="position in the ECEF frame",
                                units=['m', 'm', 'm', 'm^2', 'm^2', 'm^2', 'm^2', 'm^2', 'm^2'],
                                legend=['pos_x', 'pos_y', 'pos_z', 'cov_xx', 'cov_yy', 'cov_zz', 'cov_xy', 'cov_xz',
                                        'cov_yz'],
                                title="Estimated Position in the ECEF Frame",
                                time_cols=(0, 1),
                                data_cols=(2, 3, 4, 5, 6, 7, 8, 9, 10))

        # clock bias
        self.clock_bias = CSVData(name="clock_bias",
                                  description="Receiver Clock Bias",
                                  units=['s', 's^2'],
                                  legend=['clock_bias', 'cov'],
                                  title="Estimated Receiver Clock Bias",
                                  time_cols=(0, 1),
                                  data_cols=(2, 3))

        # clock bias rate
        self.clock_bias_rate = CSVData(name="clock_bias_rate",
                                       description="Receiver Clock Bias Rate",
                                       units=['s/s'],
                                       legend=['constellation', 'clock rate', 'cov'],
                                       title="Receiver Clock Bias Rate",
                                       time_cols=(0, 1),
                                       data_cols=(2, 3, 4))

        # velocity in ECEF
        self.velocity = CSVData(name="velocity",
                                description="velocity in the ECEF frame",
                                units=['m/s', 'm/s', 'm/s'],
                                legend=['vel_x', 'vel_y', 'vel_z'],
                                title="Velocity (ECEF)",
                                time_cols=(0, 1),
                                data_cols=(2, 3, 4))

        # Dilution of Precision (DOP)
        self.dop_ecef = CSVData(name="DOP_ECEF",
                                description="Dilution of Precision in ECEF frame",
                                units=['', '', '', '', '', ''],
                                legend=['dop_x', 'dop_y', 'dop_z', 'dop_time', 'dop_geometry', 'dop_position'],
                                title="Dilution of Precision (ECEF)",
                                time_cols=(0, 1),
                                data_cols=(2, 3, 4, 5, 6, 7))

        self.dop_local = CSVData(name="DOP_local",
                                 description="Dilution of Precision in local (ENU) frame",
                                 units=['', '', '', ''],
                                 legend=['dop_east', 'dop_north', 'dop_up', 'dop_horizontal'],
                                 title="Dilution of Precision (ENU)",
                                 time_cols=(0, 1),
                                 data_cols=(2, 3, 4, 5))

        # prefit residuals
        self.prefit_residuals = CSVData(name="prefit_residuals",
                                        description="Pseudorange Prefit Residuals",
                                        title="Pseudorange Prefit Residuals",
                                        time_cols=(0, 1),
                                        data_cols=(2, 3, 4, 5))

        # postfit residuals
        self.postfit_residuals = CSVData(name="postfit_residuals",
                                         description="Pseudorange Postfit Residuals",
                                         title="Pseudorange Postfit Residuals",
                                         time_cols=(0, 1),
                                         data_cols=(2, 3, 4, 5))

        # velocity prefit residuals
        self.vel_prefit_residuals = CSVData(name="vel_prefit_residuals",
                                            description="Pseudorange Rate Prefit Residuals",
                                            units=None,
                                            legend=None,
                                            title="Pseudorange Prefit Residuals",
                                            time_cols=(0, 1),
                                            data_cols=(2, 3, 4, 5))

        # postfit residuals
        self.vel_postfit_residuals = CSVData(name="vel_postfit_residuals",
                                             description="Pseudorange Rate Postfit Residuals",
                                             units=None,
                                             legend=None,
                                             title="Pseudorange Rate Postfit Residuals",
                                             time_cols=(0, 1),
                                             data_cols=(2, 3, 4, 5))

        # satellite azimuth elevation
        self.satellite_azel = CSVData(name="satellite_azel",
                                      description="Satellite Azimuth Elevation",
                                      units=['deg', 'deg'],
                                      legend=['sat', 'azimuth', 'elevation'],
                                      title="Satellite Azimuth and Elevation",
                                      time_cols=(0, 1),
                                      data_cols=(2, 3, 4),
                                      )

        # Inter System Bias ISB
        self.isb = CSVData(name="isb",
                           description="Inter System Bias ISB",
                           units=['s'],
                           title="Inter System Bias ISB",
                           time_cols=(0, 1),
                           data_cols=(2, 3),
                           )

        # Iono
        self.iono = CSVData(name="iono",
                            description="Ionosphere",
                            units=['m'],
                            legend=['sat', 'iono', 'cov'],
                            title="Estimated Ionosphere Delay",
                            time_cols=(0, 1),
                            data_cols=(2, 3, 4),
                            )

        self.tropo_wet = CSVData(name="tropo",
                                 description="Troposphere",
                                 units=['m'],
                                 legend=['tropo_wet'],
                                 title="Estimated Troposphere (Wet Delay)",
                                 time_cols=(0, 1),
                                 data_cols=(2, 3),
                                 )
        # available data for the current simulation
        self._available = []

    def __str__(self):
        return f'{type(self).__name__}( DataManager for GNSS Run Performance Evaluation )'

    def get_data(self, name):
        """
        Returns the required CSVData object from the data manager.

        Arguments:
            name(str): name of CSVData object to fetch

        Returns:
            CSVData: returns the fetched CSVData object.
                If the requested data is unavailable, raises a ValueError exception
        """
        if name in self._available:
            return getattr(self, name)
        else:
            raise ValueError(f'{name} is not available.')

    @property
    def available(self):
        return self._available

    def read_data(self, output_folder):
        """ TODO: add docstrings"""
        for output_name in GnssRunStorageManager.__rx_based__ + GnssRunStorageManager.__sat_based__:
            try:
                print("Reading", output_name)
                getattr(self, output_name).read_data(output_folder / OUTPUT_FILENAME_MAP[output_name])
                self._available.append(output_name)
            except Exception as e:
                print(f"[WARN] did not successfully import data {output_name} due to: {e}")
