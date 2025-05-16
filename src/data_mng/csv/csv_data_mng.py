""" Module with CSV Data Manager that manages all input data for the Post Processing GNSS Run """
import logging

from src.io.states import OUTPUT_FILENAME_MAP
from src.data_mng.container import Container

from .csv_data import CSVData


class GnssRunStorageManager(Container):
    """
        Data Manager class that stores all the data saved during the execution of the PNT GNSS algorithm. This
        data manager is used by the Post-Processing / Performance Evaluation algorithm.

        Derives from the :py:class:`Container` class.

        All outputs are stored as instances of the :py:class:`CSVData` class.

        Attributes:
            time(CSVData)
            position(CSVData)
            velocity(CSVData)
            clock_bias(CSVData)
            dop_ecef(CSVData)
            dop_local(CSVData)
            clock_bias_rate(CSVData)
            isb(CSVData)
            tropo_wet(CSVData)
            pr_prefit_residuals(CSVData)
            pr_postfit_residuals(CSVData)
            cp_prefit_residuals(CSVData)
            cp_postfit_residuals(CSVData)
            iono(CSVData)
            satellite_azel(CSVData)
            pr_rate_prefit_residuals(CSVData)
            pr_rate_postfit_residuals(CSVData)
            obs(CSVData)
            ambiguity(CSVData)
            phase_bias(CSVData)
            mw_obs(CSVData)
            gf_obs(CSVData)
        """
    __inputs__ = ["time", "position", "velocity", "clock_bias", "dop_ecef", "dop_local", "clock_bias_rate", "isb",
                  "tropo_wet", "pr_prefit_residuals", "pr_postfit_residuals", "iono", "satellite_azel",
                  "pr_rate_prefit_residuals", "pr_rate_postfit_residuals", "obs", "ambiguity", "phase_bias",
                  "cp_prefit_residuals", "cp_postfit_residuals", "mw_obs", "gf_obs"]
    __slots__ = __inputs__ + ["_available", "log"]

    def __init__(self, log: logging.Logger):
        super().__init__()
        self.log = log

        # time
        self.time = CSVData(name="time",
                            description="Time in GPS Time Scale",
                            title="GPS Time Scale",
                            time_cols=(0, 1),
                            data_cols=(2,))

        # position in ECEF
        self.position = CSVData(name="position",
                                description="Position in the ECEF frame",
                                units=['m', 'm', 'm', 'm^2', 'm^2', 'm^2', 'm^2', 'm^2', 'm^2'],
                                legend=['x', 'y', 'z', 'cov_xx', 'cov_yy', 'cov_zz', 'cov_xy', 'cov_xz',
                                        'cov_yz'],
                                title="Estimated Position in the ECEF Frame",
                                time_cols=(0, 1),
                                data_cols=(2, 3, 4, 5, 6, 7, 8, 9, 10))

        # clock bias
        self.clock_bias = CSVData(name="clock_bias",
                                  description="Receiver Clock Bias",
                                  units=['s', 's^2'],
                                  # legend=['clock_bias', 'cov'],
                                  title="Estimated Receiver Clock Bias",
                                  time_cols=(0, 1),
                                  data_cols=(2, 3))

        # clock bias rate
        self.clock_bias_rate = CSVData(name="clock_bias_rate",
                                       description="Receiver Clock Bias Rate",
                                       units=['', 's/s', '(s/s)^2'],
                                       legend=['constellation', 'clock rate', 'cov'],
                                       title="Receiver Clock Bias Rate",
                                       time_cols=(0, 1),
                                       data_cols=(2, 3, 4))

        # velocity in ECEF
        self.velocity = CSVData(name="velocity",
                                description="Velocity in the ECEF frame",
                                units=['m/s', 'm/s', 'm/s', '(m/s)^2', '(m/s)^2', '(m/s)^2', '(m/s)^2', '(m/s)^2',
                                       '(m/s)^2'],
                                legend=['x', 'y', 'z', 'cov_xx', 'cov_yy', 'cov_zz', 'cov_xy', 'cov_xz',
                                        'cov_yz'],
                                title="Estimated Velocity in the ECEF Frame",
                                time_cols=(0, 1),
                                data_cols=(2, 3, 4, 5, 6, 7, 8, 9, 10))

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

        # pseudorange prefit residuals
        self.pr_prefit_residuals = CSVData(name="pr_prefit_residuals",
                                           description="Pseudorange Prefit Residuals",
                                           title="Pseudorange Prefit Residuals",
                                           time_cols=(0, 1),
                                           data_cols=(2, 3, 4, 5),
                                           func_filter=lambda df: df[df["data_type"].str.contains("PR")])

        # pseudorange postfit residuals
        self.pr_postfit_residuals = CSVData(name="pr_postfit_residuals",
                                            description="Pseudorange Postfit Residuals",
                                            title="Pseudorange Postfit Residuals",
                                            time_cols=(0, 1),
                                            data_cols=(2, 3, 4, 5),
                                            func_filter=lambda df: df[df["data_type"].str.contains("PR")])

        # carrier phase prefit residuals
        self.cp_prefit_residuals = CSVData(name="cp_prefit_residuals",
                                           description="Carrier Phase Prefit Residuals",
                                           title="Carrier Phase Prefit Residuals",
                                           time_cols=(0, 1),
                                           data_cols=(2, 3, 4, 5),
                                           func_filter=lambda df: df[df["data_type"].str.contains("CP")])

        # carrier phase postfit residuals
        self.cp_postfit_residuals = CSVData(name="cp_postfit_residuals",
                                            description="Carrier Phase Postfit Residuals",
                                            title="Carrier Phase Postfit Residuals",
                                            time_cols=(0, 1),
                                            data_cols=(2, 3, 4, 5),
                                            func_filter=lambda df: df[df["data_type"].str.contains("CP")])

        # pr rate prefit residuals
        self.pr_rate_prefit_residuals = CSVData(name="pr_rate_prefit_residuals",
                                                description="Pseudorange Rate Prefit Residuals",
                                                title="Pseudorange Rate Prefit Residuals",
                                                time_cols=(0, 1),
                                                data_cols=(2, 3, 4, 5))

        # pr rate postfit residuals
        self.pr_rate_postfit_residuals = CSVData(name="pr_rate_postfit_residuals",
                                                 description="Pseudorange Rate Postfit Residuals",
                                                 title="Pseudorange Rate Postfit Residuals",
                                                 time_cols=(0, 1),
                                                 data_cols=(2, 3, 4, 5))

        # satellite azimuth elevation
        self.satellite_azel = CSVData(name="satellite_azel",
                                      description="Satellite Azimuth Elevation",
                                      units=['', 'deg', 'deg'],
                                      legend=['sat', 'azimuth', 'elevation'],
                                      title="Satellite Azimuth and Elevation",
                                      time_cols=(0, 1),
                                      data_cols=(2, 3, 4),
                                      )

        # Inter System Bias ISB
        self.isb = CSVData(name="isb",
                           description="Inter System Bias ISB",
                           units=['s', 's^2'],
                           title="Inter System Bias ISB",
                           time_cols=(0, 1),
                           data_cols=(2, 3),
                           )

        # Iono
        self.iono = CSVData(name="iono",
                            description="Ionosphere",
                            units=['', 'm', 'm^2'],
                            legend=['sat', 'iono', 'cov'],
                            title="Estimated Ionosphere Delay",
                            time_cols=(0, 1),
                            data_cols=(2, 3, 4),
                            )

        self.tropo_wet = CSVData(name="tropo",
                                 description="Troposphere",
                                 units=['m', 'm^2'],
                                 legend=['tropo_wet', 'cov'],
                                 title="Estimated Troposphere (Wet Delay)",
                                 time_cols=(0, 1),
                                 data_cols=(2, 3),
                                 )

        # GNSS observation data
        self.obs = CSVData(name="obs",
                           description="GNSS Observations",
                           title="GNSS Observations",
                           time_cols=(0, 1),
                           data_cols=(2, 3, 4, 5))

        # GNSS observation data (Melbourne-Wubbena Combinations)
        self.mw_obs = CSVData(name="mw_obs",
                              description="Melbourne-Wubbena Observations",
                              title="Melbourne-Wubbena Observations",
                              time_cols=(0, 1),
                              data_cols=(2, 3, 4, 5))

        # GNSS observation data (Geometry-Free Combinations)
        self.gf_obs = CSVData(name="gf_obs",
                              description="Geometry-Free Observations",
                              title="Geometry-Free Observations",
                              time_cols=(0, 1),
                              data_cols=(2, 3, 4, 5))

        # Ambiguity
        self.ambiguity = CSVData(name="ambiguity",
                                 description="Ambiguity",
                                 units=['', '', 'm', 'm^2'],
                                 legend=['sat', 'datatype', 'ambiguity', 'cov'],
                                 title="Estimated Ambiguity",
                                 time_cols=(0, 1),
                                 data_cols=(2, 3, 4, 5),
                                 )

        # Phase Bias
        self.phase_bias = CSVData(name="phase_bias",
                                  description="Receiver Phase Bias",
                                  units=['', '', 'm', 'm^2'],
                                  legend=['constellation', 'datatype', 'phase_bias', 'cov'],
                                  title="Estimated Receiver Phase Bias",
                                  time_cols=(0, 1),
                                  data_cols=(2, 3, 4, 5),
                                  )

        # available data for the current simulation
        self._available = []

    def __str__(self):
        return f'{type(self).__name__}( DataManager for GNSS Run Performance Evaluation )'

    def get_data(self, name):
        """
        Returns the required CSVData object from the data manager.

        Args:
            name(str): name of CSVData object to fetch

        Returns:
            CSVData: returns the fetched CSVData object.

        Raises:
            ValueError: If the requested data is unavailable, an exception is raised
        """
        if name in self._available:
            return getattr(self, name)
        else:
            raise ValueError(f'{name} is not available.')

    @property
    def available(self):
        return self._available

    def read_data(self, output_folder):
        """ Read all input files available for this run """
        for output_name in GnssRunStorageManager.__inputs__:
            try:
                in_file = output_folder / OUTPUT_FILENAME_MAP[output_name]
                self.log.info(f"reading data '{output_name}' stored in file {in_file}...")

                getattr(self, output_name).read_data(in_file)
                self._available.append(output_name)
            except Exception as e:
                self.log.warning(f"Did not successfully import data '{output_name}' due to: {repr(e)}")
