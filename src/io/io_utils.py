"""" Module with common utilities for Input / Output file operations"""
import time
import os
import re
from datetime import datetime

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


def get_last_run_folder():
    """
    Get the most recent run folder in the RUNS_PATH directory.

    Returns:
        str: the name of the most recent run folder
    """

    home_dir = os.path.expanduser(RUNS_PATH)

    # Regular expression to match the timestamp format YYYY-MM-DDTHHMMSS
    time_pattern = re.compile(r"(\d{4}-\d{2}-\d{2}T\d{2}H\d{2}M\d{2}S)")

    # Get all valid run folders
    run_folders = [
        f for f in os.listdir(home_dir)
        if os.path.isdir(os.path.join(home_dir, f)) and time_pattern.fullmatch(f)
    ]

    # Sort folders by timestamp (latest first)
    run_folders.sort(key=lambda x: datetime.strptime(x, "%Y-%m-%dT%HH%MM%SS"), reverse=True)

    # Get the last (most recent) run folder
    last_run_folder = run_folders[0] if run_folders else None

    return last_run_folder
