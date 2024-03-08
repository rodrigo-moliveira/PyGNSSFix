# from src.constants import OUTPUT_FILENAME_MAP
# from src.data_mng.container import Container
# from src.io.csv_data import CSVData, read_csv_file
#
#
# # TODO: quando chegar ao INS, criar um class mae com as features generices e depois separar em duas
#
# def parse_satellite_data_file(data_in, cols):
#     data_out = {}
#     for line in data_in:
#
#         # fetch data for this line
#         epoch = (line[0], line[1])  # week and sow
#         sat = line[2]  # satellite ID
#         data = []
#         for col in cols:
#             data.append(float(line[col]))
#
#         # add epoch to outer dict
#         if epoch not in data_out.keys():
#             data_out[epoch] = {}
#
#         # add this data point
#         data_out[epoch][sat] = data
#
#     return data_out
#
#
# class GnssRunStorageManager(Container):
#     __slots__ = ["time", "position", "velocity", "prefit_residuals", "clock_bias",
#                  "postfit_residuals", "dop_ecef", "dop_local", "satellite_azel", "_available"]
#
#     def __init__(self):
#         super().__init__()
#
#         # time
#         self.time = CSVData(name="time",
#                             description="Time steps",
#                             units=[''],
#                             legend=[''],
#                             title="Time")
#
#         # position in ECEF
#         self.position = CSVData(name="position",
#                                 description="position in the ECEF frame",
#                                 units=['m', 'm', 'm'],
#                                 legend=['pos_x', 'pos_y', 'pos_z'],
#                                 title="Position (ECEF)")
#
#         # clock bias
#         self.clock_bias = CSVData(name="clock_bias",
#                                   description="Receiver Clock Bias",
#                                   units=['s'],
#                                   legend=['time'],
#                                   title="Receiver Clock Bias")
#
#         # velocity in ECEF
#         self.velocity = CSVData(name="velocity",
#                                 description="velocity in the ECEF frame",
#                                 units=['m/s', 'm/s', 'm/s'],
#                                 legend=['vel_x', 'vel_y', 'vel_z'],
#                                 title="Velocity (ECEF)")
#
#         # Dilution of Precision (DOP)
#         self.dop_ecef = CSVData(name="DOP_ECEF",
#                                 description="Dilution of Precision in ECEF frame",
#                                 units=['', '', '', ''],
#                                 legend=['', '', '', ''],
#                                 title="Dilution of Precision (ECEF)")
#
#         self.dop_local = CSVData(name="DOP_local",
#                                  description="Dilution of Precision in local (ENU) frame",
#                                  units=['', '', '', ''],
#                                  legend=['', '', '', ''],
#                                  title="Dilution of Precision (ENU)")
#
#         # prefit residuals
#         self.prefit_residuals = CSVData(name="prefit_residuals",
#                                         description="Prefit Residuals",
#                                         units=['m'],
#                                         legend=[' '],
#                                         title="Prefit Residuals")
#
#         # postfit residuals
#         self.postfit_residuals = CSVData(name="postfit_residuals",
#                                          description="Postfit Residuals",
#                                          units=['m'],
#                                          legend=[' '],
#                                          title="Postfit Residuals")
#
#         # satellite azimuth elevation
#         self.satellite_azel = CSVData(name="satellite_azel",
#                                       description="Satellite Azimuth Elevation",
#                                       units=['deg', 'deg'],
#                                       legend=['azimuth', 'elevation'],
#                                       title="Satellite Azimuth and Elevation")
#
#         # available data for the current simulation
#         self._available = []
#
#     def __str__(self):
#         return f'{type(self).__name__}( DataManager for GNSS Run Performance Evaluation )'
#
#     def add_data(self, data_name, data, units=None):
#         """
#         Add data to available.
#         Args:
#             data_name: data name, str
#             data: a scalar, a numpy array or a dict of the above two. If data is a dict, each
#                 value in it should be of same type (scalar or numpy array), same size and same
#                 units.
#             units: Units of the data. If you know clearly no units convertion is needed, set
#                 units to None. If you do not know what units are used in the class InsDataMgr,
#                 you'd better provide the units of the data. Units convertion will be done
#                 automatically.
#                 If data is a scalar, units should be a list of one string to define its unit.
#                 If data is a numpy of size(m,n), units should be a list of n strings
#                 to define the units.
#         """
#         if data_name in self.__slots__:
#             sim = getattr(self, data_name, None)
#             if sim is not None:
#                 sim.add_data(data, units)
#
#                 # add to 'available' list
#                 if data_name not in self._available:
#                     self._available.append(data_name)
#         else:
#             raise ValueError(f"Unsupported data: {data_name}, not in {self.__slots__}")
#
#     def get_data(self, data):
#         """
#         Get data section of data_names.
#         Args:
#             data: name of data to get
#         Returns:
#             data: a list of data corresponding to data_names.
#             If there is any unavailable data in data_names, return None
#         """
#         if data in self._available:
#             return getattr(self, data)
#         else:
#             raise ValueError(f'{data} is not available.')
#
#     @property
#     def available(self):
#         return self._available
#
#     def read_data(self, output_folder):
#         time = read_csv_file(output_folder / OUTPUT_FILENAME_MAP["time"],
#                              header=True, columns=(0, 1))
#         position = read_csv_file(output_folder / OUTPUT_FILENAME_MAP["position"],
#                                  header=True, columns=(2, 3, 4))
#         dop_ecef = read_csv_file(output_folder / OUTPUT_FILENAME_MAP["dop_ecef"],
#                                  header=True, columns=(2, 3, 4, 5, 6, 7))
#         dop_enu = read_csv_file(output_folder / OUTPUT_FILENAME_MAP["dop_local"],
#                                 header=True, columns=(2, 3, 4, 5))
#         clock_bias = read_csv_file(output_folder / OUTPUT_FILENAME_MAP["clock_bias"],
#                                    header=True, columns=(2,))
#         satellite_azel = read_csv_file(output_folder / OUTPUT_FILENAME_MAP["satellite_azel"],
#                                        header=True, columns=(0, 1, 2, 3, 4), dtype=str)
#         residuals = read_csv_file(output_folder / OUTPUT_FILENAME_MAP["postfit_residuals"],
#                                   header=True, columns=(0, 1, 2, 3), dtype=str)
#
#         self.add_data("time", time)
#         self.add_data("position", position)
#         self.add_data("dop_ecef", dop_ecef)
#         self.add_data("dop_local", dop_enu)
#         self.add_data("clock_bias", clock_bias)
#         self.add_data("satellite_azel", parse_satellite_data_file(satellite_azel, cols=(3, 4)))
#         self.add_data("postfit_residuals", parse_satellite_data_file(residuals, cols=(3,)))
