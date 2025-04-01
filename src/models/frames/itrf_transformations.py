"""Module to perform transformations between different realizations of the International Terrestrial Reference
Frame (ITRF).
"""
import json
import jsonschema
from pyproj import Transformer, CRS, database
import numpy as np

from src import PROJECT_PATH
from src.errors import ITRFError
from src.data_types.date import Epoch


class ITRF_Transformation:
    """
    Class to perform transformations between different realizations of the International Terrestrial Reference
    Frame (ITRF).

    The transformations are based on Helmert transformation parameters taken from the EPSG Geodetic
    Parameter Dataset.

    Two methods are available to perform the transformation:
    1. Using the pyproj library to perform the transformation directly.
    2. Using Helmert transformation parameters from a JSON file (provided by the user) to perform the transformation.

    For more information, consult https://epsg.io/

    Attributes:
        in_frame (str): Name of the input frame.
        out_frame (str): Name of the output frame.
        transformer (pyproj.Transformer): Pyproj transformer object.
        json_data (dict): Dictionary containing transformation parameters from a JSON
        forward (bool): Flag to indicate the direction of the transformation (default: True).
    """
    _cache = {}

    def __init__(self, in_frame: str, out_frame: str):
        """
        Constructor for the ITRF_Transformation class.

        Args:
            in_frame (str): Name of the input frame
            out_frame (str): Name of the output frame.
        """
        self.in_frame = in_frame
        self.out_frame = out_frame
        self.transformer = None
        self.json_data = None
        self.forward = True

    @staticmethod
    def factory(in_frame: str, out_frame: str, file_path=None):
        """
        Factory method to create an `ITRF_Transformation` object.

        If the transformation parameters (json file) are not provided, the method will use the pyproj library to
        fetch the transformation parameters from the EPSG Geodetic Parameter Dataset.

        Args:
            in_frame (str): Name of the input frame.
            out_frame (str): Name of the output frame.
            file_path (str): Path to the JSON file containing transformation parameters (optional).

        Returns:
            ITRF_Transformation: An `ITRF_Transformation` object.
        """
        try:
            if f"{in_frame}_{out_frame}" not in ITRF_Transformation._cache:
                from src.common_log import get_logger, IO_LOG
                log = get_logger(IO_LOG)

                obj = ITRF_Transformation(in_frame, out_frame)
                if file_path is None:
                    log.info("Constructing ITRF Transformation class. "
                             "Using pyproj to fetch transformation parameters...")
                    obj.from_crs(log)
                else:
                    log.info("Constructing ITRF Transformation class. "
                             f"Using JSON file {file_path} to fetch transformation parameters...")
                    obj.from_json(file_path, log)
                ITRF_Transformation._cache[f"{in_frame}_{out_frame}"] = obj
        except Exception as e:
            raise ITRFError(f"Error in ITRF Transformation factory: {e}")

        return ITRF_Transformation._cache[f"{in_frame}_{out_frame}"]

    def from_crs(self, log):
        """
        Fetches the transformation parameters from the EPSG Geodetic Parameter Dataset using the pyproj library.

        Initializes the `transformer` attribute with the transformation parameters.

        Args:
            log (logging.Logger): Logger object.
        """

        # Query the EPSG Geodetic Parameter Dataset
        search_results = database.query_crs_info(auth_name="EPSG")

        # Filter results for the input and output frames
        search_in_crs = [crs for crs in search_results if self.in_frame == crs.name]
        if len(search_in_crs) == 0:
            raise ValueError(f"Input frame {self.in_frame} not found in EPSG database.")
        in_frame = search_in_crs[0]

        search_out_crs = [crs for crs in search_results if self.out_frame == crs.name]
        if len(search_out_crs) == 0:
            raise ValueError(f"Output frame {self.out_frame} not found in EPSG database.")
        out_frame = search_out_crs[0]

        # Log the input and output frames
        log.info(f"[Input frame for ITRF Transformations] EPSG Code: {in_frame.code}, Name: {in_frame.name}")
        log.info(f"[Output frame for ITRF Transformations] EPSG Code: {out_frame.code}, Name: {out_frame.name}")

        crs_from = CRS.from_epsg(in_frame.code)
        crs_to = CRS.from_epsg(out_frame.code)

        # Get available transformations
        self.transformer = Transformer.from_crs(crs_from, crs_to, always_xy=True)

    def from_json(self, file_path, log):
        """
        Loads the transformation parameters from a JSON file.

        Validates the JSON file against the schema and checks if the input and output frames match the user-defined
        frames.
        Initializes the `json_data` attribute with the transformation parameters.

        Args:
            file_path (str): Path to the JSON file.
            log (logging.Logger): Logger object.
        """
        # Load the JSON file and schema
        schema_file = PROJECT_PATH / "src/io/config/resources/itrf_epsg.schema.json"
        json_data = ITRF_Transformation.load_json_file(file_path)
        json_schema = ITRF_Transformation.load_json_file(schema_file)

        # Validate the JSON
        jsonschema.validate(instance=json_data, schema=json_schema)

        if json_data["source_crs"]["name"] == self.in_frame and json_data["target_crs"]["name"] == self.out_frame:
            pass
        elif json_data["source_crs"]["name"] == self.out_frame and json_data["target_crs"]["name"] == self.in_frame:
            # backward transformation
            self.forward = False
        else:
            raise ValueError(f"User-defined input frame {self.in_frame} and output frame {self.out_frame} "
                             f"do not match JSON source CRS ({json_data['source_crs']['name']}) and "
                             f"target CRS ({json_data['target_crs']['name']}).")

        log.info(f"Successfully loaded validated JSON file {file_path} for ITRF transformation.")
        log.info(f"[Input frame for ITRF Transformations] Name: {json_data['source_crs']['name']}")
        log.info(f"[Output frame for ITRF Transformations] Name: {json_data['target_crs']['name']}")
        self.json_data = json_data

    def transform(self, input_position, time: float or Epoch):
        """
        Applies the transformation to the input coordinates.

        Args:
            input_position (numpy.ndarray or list): Input coordinates in the input frame.
            time (float or Epoch): Target epoch, either in decimal years (floating point) or as an Epoch object.

        Returns:
            numpy.ndarray: Transformed coordinates in the output frame.
        """
        if len(input_position) != 3:
            raise ValueError("Input position must be a 3-element list or numpy array.")

        if isinstance(time, Epoch):
            time = time.decimal_year

        if self.transformer is not None:
            result = self.transformer.transform(xx=input_position[0], yy=input_position[1], zz=input_position[2],
                                                tt=time)
            position = np.array(result[0:3])

        elif self.json_data is not None:
            # Get transformation parameters from JSON
            T, R_matrix, S = ITRF_Transformation.get_transformation_parameters(time, self.json_data)

            # Apply Helmert transformation
            position = self.helmert_transform(input_position, T, R_matrix, S)

        else:
            raise ITRFError("No transformation data available.")
        return position

    @staticmethod
    def get_transformation_parameters(epoch, data):
        """
        Computes transformation parameters for a given epoch.

        Args:
            epoch (float): Target epoch (decimal year).
            data (dict): Dictionary parsed from JSON file containing transformation parameters.

        Returns:
            tuple: Helmert transformation parameters (T, R_matrix, S).
        """
        # Read parameters from JSON
        epoch_ref, T_ref, R_ref, S_ref, T_rate, R_rate, S_rate = ITRF_Transformation.read_transformation_json(data)

        # Compute epoch difference
        factor = ITRF_Transformation._extract_conversion_factor(data, "Parameter reference epoch")
        dt = epoch * factor - epoch_ref

        # Compute adjusted parameters
        T = T_ref + T_rate * dt  # Adjusted translation (meters)
        R_rad = R_ref + R_rate * dt  # Adjusted rotation (radians)
        S = S_ref + S_rate * dt  # Adjusted scale factor (unit-less)

        # Compute rotation matrix using small-angle approximation
        R_matrix = np.array([
            [1, -R_rad[2], R_rad[1]],
            [R_rad[2], 1, -R_rad[0]],
            [-R_rad[1], R_rad[0], 1]
        ])

        return T, R_matrix, S

    def helmert_transform(self, X, T, R, S):
        """
        Applies the Helmert transformation.

        Args:
            X (numpy.ndarray): (3,) Original coordinates in meters.
            T (numpy.ndarray): (3,) Translation vector in meters.
            R (numpy.ndarray): (3,3) Rotation matrix.
            S (float): Scale factor (unit-less, converted from parts per billion).

        Returns:
            numpy.ndarray: Transformed coordinates in meters.
        """
        if self.forward:
            return T + (1 + S) * (R @ X)
        else:
            return 1/(1 + S) * (R.T @ (X - T))

    @staticmethod
    def _extract_param(data, name):
        """ Helper function to extract a parameter's value and convert using its unit factor."""
        param = next(param for param in data["parameters"] if param["name"] == name)
        return param["value"] * param["unit"]["conversion_factor"]

    @staticmethod
    def _extract_conversion_factor(data, name):
        """ Helper function to extract a parameter's conversion factor."""
        param = next(param for param in data["parameters"] if param["name"] == name)
        return param["unit"]["conversion_factor"]

    @staticmethod
    def read_transformation_json(data):
        """
        Reads the transformation parameters from the provided JSON structure.

        Args:
            data (dict): Dictionary parsed from a JSON file containing transformation parameters.

        Returns:
            tuple: (epoch_ref, T_ref, R_ref, S_ref, T_rate, R_rate, S_rate)
        """

        # Extract reference epoch
        epoch_ref = ITRF_Transformation._extract_param(data, "Parameter reference epoch")

        # Extract translation (converted from mm to meters)
        T_ref = np.array(
            [ITRF_Transformation._extract_param(data, f"{axis}-axis translation") for axis in ["X", "Y", "Z"]])

        # Extract rotation (converted from milli arc-seconds to radians)
        R_ref = np.array([ITRF_Transformation._extract_param(data, f"{axis}-axis rotation")
                          for axis in ["X", "Y", "Z"]])

        # Extract scale difference (converted from ppb to unit-less)
        S_ref = ITRF_Transformation._extract_param(data, "Scale difference")

        # Extract rates of change
        T_rate = np.array(
            [ITRF_Transformation._extract_param(data, f"Rate of change of {axis}-axis translation") for axis in
             ["X", "Y", "Z"]])
        R_rate = np.array(
            [ITRF_Transformation._extract_param(data, f"Rate of change of {axis}-axis rotation") for axis in
             ["X", "Y", "Z"]])
        S_rate = ITRF_Transformation._extract_param(data, "Rate of change of Scale difference")

        return epoch_ref, T_ref, R_ref, S_ref, T_rate, R_rate, S_rate

    @staticmethod
    def load_json_file(file_path):
        with open(file_path, "r", encoding="utf-8") as file:
            return json.load(file)


# Example usage
if __name__ == "__main__":
    # Define input coordinates (X, Y, Z) in meters and time in decimal years
    _X = 4027894.0060
    _Y = 307045.6000
    _Z = 4919474.9100
    _time = 2000.00

    # Unit test 1) Forward Transformation
    tf_crs = ITRF_Transformation.factory("ITRF93", "ITRF2020")
    result_crs = tf_crs.transform([_X, _Y, _Z], _time)
    ITRF_Transformation._cache.clear()
    tf_json = ITRF_Transformation.factory("ITRF93", "ITRF2020",
                                          file_path="../../../workspace/geo_time_data/ITRF/EPSG_9998.json")
    result_json = tf_json.transform([_X, _Y, _Z], _time)

    # Report
    print("Unit Test 1: Forward Transformation")
    print(f"Transformed position from CRS: {result_crs}")
    print(f"Transformed position from JSON: {result_json}")
    print(f"Transformed position from online tool: [4027894.0539, 307045.5594, 4919474.9073]\n")

    # Unit test 2) Backward Transformation
    tf_crs = ITRF_Transformation.factory("ITRF2020", "ITRF93")
    result_crs = tf_crs.transform([_X, _Y, _Z], _time)
    ITRF_Transformation._cache.clear()
    tf_json = ITRF_Transformation.factory("ITRF2020", "ITRF93",
                                          file_path="../../../workspace/geo_time_data/ITRF/EPSG_9998.json")
    result_json = tf_json.transform([_X, _Y, _Z], _time)

    # Report
    print("Unit Test 2: Backward Transformation")
    print(f"Transformed position from CRS: {result_crs}")
    print(f"Transformed position from JSON: {result_json}")


"""
Validation Example from https://epncb.oma.be/_productsservices/coord_trans/

#_Station Frame    Epoch        X[m]          Y[m]          Z[m]      VX[m/yr] VY[m/yr] VZ[m/yr]
Station_1 ITRF93   2000.00  4027894.0060   307045.6000  4919474.9100
Station_1 ITRF2020 2000.00  4027894.0539   307045.5594  4919474.9073

Station_2 ITRF93   2000.00  4027894.0060   307045.6000  4919474.9100   0.0100   0.2000   0.0300
Station_2 ITRF2020 2000.00  4027894.0539   307045.5594  4919474.9073   0.0170   0.1962   0.02

Station_3 ITRF93   2025.75  4027894.0060   307045.6000  4919474.9100
Station_3 ITRF2020 2025.75  4027894.2329   307045.4608  4919474.8600

Station_4 ITRF93   2025.75  4027894.0060   307045.6000  4919474.9100   0.0100   0.2000   0.0300
Station_4 ITRF2020 2025.75  4027894.2329   307045.4608  4919474.8600   0.0170   0.1962   0.0282
"""
