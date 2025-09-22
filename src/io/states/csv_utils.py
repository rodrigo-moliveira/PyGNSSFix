# """ Module with utility functions to read CSV files """
# import numpy as np
#
#
# def parse_satellite_data_file(data_in, cols):
#     """
#     Parses satellite data from an input source and organizes it into a nested dictionary.
#
#     Parameters:
#         data_in(list[list]): A list of lists where each sublist represents a line of satellite data.
#             Each sublist is expected to contain:
#             - week (int): The GPS week number.
#             - sow (float): The seconds of the week.
#             - sat (str): The satellite ID.
#             - Other values: The satellite data for various parameters.
#         cols(list[int] or None): A list of column indices specifying which columns of data_in contain the desired
#             satellite data to be parsed and included in the output.
#
#     Returns:
#         dict: A nested dictionary where the outer keys are tuples (epoch) containing the week and sow, and the inner
#             keys are satellite IDs. The values are lists of floats representing the satellite data corresponding to
#             the specified columns.
#
#     Examples:
#     ---------
#     Parse satellite data with columns of interest:
#
#     >>> data_in_example = [
#     ...     [1920, 432000, 'G01', 0.1, 0.2, 0.3, 0.4],
#     ...     [1920, 432000, 'G02', 0.5, 0.6, 0.7, 0.8]
#     ... ]
#     >>> cols_example = [3, 4, 5]
#     >>> parsed_data = parse_satellite_data_file(data_in_example, cols_example)
#     >>> print(parsed_data)
#     {
#         (1920, 432000): {
#             'G01': [0.1, 0.2, 0.3],
#             'G02': [0.5, 0.6, 0.7]
#         }
#     }
#     """
#     print(data_in)
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
#     print(data_out)
#     exit()
#     return data_out
#
#
# def read_csv_file(filepath, **options):
#     """
#     Reads a CSV file and returns its contents as a NumPy array.
#
#     Parameters:
#         filepath(pathlib.Path or str): path to the input CSV file
#         **options(dict): optional
#             Additional options for reading the CSV file. The dictionary may contain the following keys:
#             - 'cols' (tuple of int, optional): Specifies the column indices to be read from the CSV file.
#                     If not provided, all columns are read.
#             - 'header' (bool, optional): Indicates whether the first line of the CSV file should be skipped
#                     (assumed to be the header). Default is True.
#             - 'dtype' (data-type, optional): The desired data type of the returned array. Default is float.
#
#     Returns:
#         numpy.ndarray: A NumPy array containing the data read from the CSV file.
#     """
#     cols = options.get("cols", None)
#     skip_first = 1 if options.get("header", True) else 0
#     dtype = options.get("dtype", float)
#
#     data = np.loadtxt(filepath,
#                       skiprows=skip_first,
#                       delimiter=',',
#                       usecols=cols,
#                       dtype=dtype)
#     return data
