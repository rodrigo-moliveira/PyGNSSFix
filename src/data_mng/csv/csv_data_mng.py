from src.io.states import OUTPUT_FILENAME_MAP
from src.data_mng.container import Container

from .csv_data import CSVData


class GnssRunStorageManager(Container):
    __rx_based__ = ["time", "position", "velocity", "clock_bias",
                    "dop_ecef", "dop_local"]
    __sat_based__ = ["prefit_residuals", "postfit_residuals"]

    __slots__ = __rx_based__ + __sat_based__ + ["_available"]

    def __init__(self):
        super().__init__()
        # TODO: colocar todos os possiveis ficheiros (incluindo ISB, Iono...)
        # time
        self.time = CSVData(name="time",
                            description="Time in GPS Time Scale",
                            units=['WN', 's'],
                            legend=['week_number', 'seconds_of_week'],
                            title="GPS Time Scale",
                            cols=(0, 1))

        # position in ECEF
        self.position = CSVData(name="position",
                                description="position in the ECEF frame",
                                units=['m', 'm', 'm'],
                                legend=['pos_x', 'pos_y', 'pos_z'],
                                title="Position (ECEF)",
                                cols=(2, 3, 4))

        # clock bias
        self.clock_bias = CSVData(name="clock_bias",
                                  description="Receiver Clock Bias",
                                  units=['s'],
                                  legend=['time'],
                                  title="Receiver Clock Bias",
                                  cols=(2,))

        # velocity in ECEF
        self.velocity = CSVData(name="velocity",
                                description="velocity in the ECEF frame",
                                units=['m/s', 'm/s', 'm/s'],
                                legend=['vel_x', 'vel_y', 'vel_z'],
                                title="Velocity (ECEF)",
                                cols=(2, 3, 4))

        # Dilution of Precision (DOP)
        self.dop_ecef = CSVData(name="DOP_ECEF",
                                description="Dilution of Precision in ECEF frame",
                                units=['', '', '', '', '', ''],
                                legend=['dop_x', 'dop_y', 'dop_z', 'dop_time', 'dop_geometry', 'dop_position'],
                                title="Dilution of Precision (ECEF)",
                                cols=(2, 3, 4, 5, 6, 7))

        self.dop_local = CSVData(name="DOP_local",
                                 description="Dilution of Precision in local (ENU) frame",
                                 units=['', '', '', ''],
                                 legend=['dop_east', 'dop_north', 'dop_up', 'dop_horizontal'],
                                 title="Dilution of Precision (ENU)",
                                 cols=(2, 3, 4, 5))

        # prefit residuals
        self.prefit_residuals = CSVData(name="prefit_residuals",
                                        description="Prefit Residuals",
                                        units=None,
                                        legend=None,
                                        title="Prefit Residuals",
                                        cols=(2, 3, 4, 5),
                                        dtype='U3')  # for residuals, we need to read the csv as string cells

        # postfit residuals
        self.postfit_residuals = CSVData(name="postfit_residuals",
                                         description="Postfit Residuals",
                                         units=None,
                                         legend=None,
                                         title="Postfit Residuals",
                                         cols=(2, 3, 4, 5),
                                         dtype='U3')  # for residuals, we need to read the csv as string cells

        # # satellite azimuth elevation
        # self.satellite_azel = CSVData(name="satellite_azel",
        #                               description="Satellite Azimuth Elevation",
        #                               units=['deg', 'deg'],
        #                               legend=['azimuth', 'elevation'],
        #                               title="Satellite Azimuth and Elevation")

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
        for output_name in GnssRunStorageManager.__rx_based__ + GnssRunStorageManager.__sat_based__:
            print("Reading", output_name)
            getattr(self, output_name).read_data(output_folder / OUTPUT_FILENAME_MAP[output_name])
            self._available.append(output_name)

        # format satellite-based datasets
        for output_name in GnssRunStorageManager.__sat_based__:
            print("Formatting satellite based", output_name)
            getattr(self, output_name).format_sat_data()
