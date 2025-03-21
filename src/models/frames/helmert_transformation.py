"""Module to perform Helmert transformation between different realizations of the International Terrestrial
 Reference Frame (ITRF).
"""

import json
import jsonschema
import numpy as np


def extract_param(data, name):
    """Helper function to extract a parameter's value and convert using its unit factor."""
    param = next(param for param in data["parameters"] if param["name"] == name)
    return param["value"] * param["unit"]["conversion_factor"]


def extract_conversion_factor(data, name):
    """Helper function to extract a parameter's conversion factor."""
    param = next(param for param in data["parameters"] if param["name"] == name)
    return param["unit"]["conversion_factor"]

def read_transformation_json(data):
    """
    Reads the transformation parameters from the provided JSON structure.

    Parameters:
        data (dict): Dictionary parsed from a JSON file containing transformation parameters.

    Returns:
        tuple: (epoch_ref, T_ref, R_ref, S_ref, T_rate, R_rate, S_rate)
    """

    # Extract reference epoch
    epoch_ref = extract_param(data, "Parameter reference epoch")

    # Extract translation (converted from mm to meters)
    T_ref = np.array([extract_param(data, f"{axis}-axis translation") for axis in ["X", "Y", "Z"]])

    # Extract rotation (converted from milliarc-seconds to radians)
    R_ref = np.array([extract_param(data, f"{axis}-axis rotation") for axis in ["X", "Y", "Z"]])

    # Extract scale difference (converted from ppb to unitless)
    S_ref = extract_param(data, "Scale difference")

    # Extract rates of change
    T_rate = np.array([extract_param(data, f"Rate of change of {axis}-axis translation") for axis in ["X", "Y", "Z"]])
    R_rate = np.array([extract_param(data, f"Rate of change of {axis}-axis rotation") for axis in ["X", "Y", "Z"]])
    S_rate = extract_param(data, "Rate of change of Scale difference")

    return epoch_ref, T_ref, R_ref, S_ref, T_rate, R_rate, S_rate


def get_transformation_parameters(epoch, data):
    """
    Computes transformation parameters for a given epoch.

    Parameters:
        epoch (float): Target epoch (decimal year).
        data (dict): Dictionary parsed from JSON file containing transformation parameters.

    Returns:
        tuple: (T, R_matrix, S)
    """
    # Read parameters from JSON
    epoch_ref, T_ref, R_ref, S_ref, T_rate, R_rate, S_rate = read_transformation_json(data)

    # Compute epoch difference
    factor = extract_conversion_factor(data, "Parameter reference epoch")
    dt = epoch*factor - epoch_ref

    # Compute adjusted parameters
    T = T_ref + T_rate * dt  # Adjusted translation (meters)
    R_rad = R_ref + R_rate * dt  # Adjusted rotation (radians)
    S = S_ref + S_rate * dt  # Adjusted scale factor (unitless)

    # Compute rotation matrix using small-angle approximation
    R_matrix = np.array([
        [1, -R_rad[2], R_rad[1]],
        [R_rad[2], 1, -R_rad[0]],
        [-R_rad[1], R_rad[0], 1]
    ])

    return T, R_matrix, S


def helmert_transform(X, T, R, S):
    """
    Applies the Helmert transformation.

    Parameters:
        X (numpy array): (3,) Original coordinates in meters (ITRF93).
        T (numpy array): (3,) Translation vector in meters.
        R (numpy array): (3,3) Rotation matrix.
        S (float): Scale factor (unitless, converted from parts per billion).

    Returns:
        numpy array: Transformed coordinates in ITRF2020.
    """
    return T + (1 + S) * (R @ X)



# Function to load a JSON file
def load_json_file(file_path):
    with open(file_path, "r", encoding="utf-8") as file:
        return json.load(file)

# Function to validate JSON against schema
def validate_json(json_data, schema):
    try:
        jsonschema.validate(instance=json_data, schema=schema)
        print("JSON is valid.")
    except jsonschema.exceptions.ValidationError as err:
        print("JSON validation error:", err)

# Example usage
if __name__ == "__main__":
    # Fetch schema from the URL
    schema_file = "../../io/config/resources/projjson.schema.json"

    # Load the JSON file
    json_data = load_json_file("../../../workspace/geo_time_data/ITRF/EPSG_9998.json")
    json_schema = load_json_file(schema_file)

    # Validate the JSON
    validate_json(json_data, json_schema)

    epoch = 2025.75  # Example epoch for transformation

    # Example Sun position in ITRF93 (from CSpice) in meters
    X_93 = np.array([4027894.0060,  307045.6000,  4919474.9100])  # Example Sun position

    # Get transformation parameters from JSON
    T, R_matrix, S = get_transformation_parameters(epoch, json_data)

    # Apply Helmert transformation
    X_2020 = helmert_transform(X_93, T, R_matrix, S)

    print("Sun position in ITRF2020:", X_2020)

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