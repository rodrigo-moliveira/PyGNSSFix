"""" Module with common utilities for Input / Output file operations"""
import time
import os

from src import RUNS_PATH


def clean_text_file(file_path):
    """
    Cleans (empties) the specified text file if it exists.

    Args:
        file_path (str): The path to the text file to be cleaned.

    Returns:
        bool: True if the file existed and was cleaned, False if the file did not exist.
    """
    try:
        with open(file_path, 'w'):
            pass  # Opening the file in 'w' mode empties the file
        return True
    except FileNotFoundError:
        return False


def create_output_dir(data_dir=None):
    """
    This function creates the output directory to save the run files.

    If the optional argument `data_dir` is passed, the function first checks if it is a valid directory/folder.
    Otherwise, the default name is used (timestamp of current time).

    Args:
        data_dir(str): directory output name (optional)
    Returns:
        str: returns the absolute path of the created directory
    Raises:
        IOError: if the creation of the output directory fails, an IOError exception is raised
    """
    # if data_dir is not specified, automatically create one
    if data_dir is None or data_dir == '' or data_dir == "default":
        data_dir = str(RUNS_PATH)
        if data_dir[-1] != '//':
            data_dir = data_dir + '//'
        data_dir = data_dir + time.strftime('%Y-%m-%dT%HH%MM%SS', time.localtime()) + '//'
        data_dir = os.path.abspath(data_dir)

    # try to create data dir
    if not os.path.exists(data_dir):
        try:
            data_dir = os.path.abspath(data_dir)
            os.makedirs(data_dir)
            os.makedirs(f"{data_dir}\\trace")
            os.makedirs(f"{data_dir}\\output")
        except:
            raise IOError(f"Cannot create dir: {data_dir}")
    return data_dir
