""" Module with useful string utility functions to support the library """

import re


def replace_whitespace_with_underscore(s):
    """
    Replace one or more whitespace characters with a single underscore.

    Args:
        s(str): input string to be processed
    Returns:
        str: the processed string, with whitespaces replaced by underscore
    """
    return re.sub(r'\s+', '_', s)


def extract_constellations(clock_string, isb_string):
    """
    Extracts the constellation names (master and slave) using regex.

    The regex patterns are the following:
        1) ``r'clock_bias\(master=(\w+)\)'``: to extract the master constellation
        2) ``r'ISB\(master=(\w+) slave=(\w+)\)'``: to extract the slave constellation

    Args:
        clock_string(str): input string in the format defined by the 1st regex expression above
        isb_string(str): input string in the format defined by the 2nd regex expression above
    Returns:
        tuple[str,str]: the extracted master and slave constellations
    """
    master, slave = "", ""
    # Define the regex pattern to match the master and slave constellations
    clock_pattern = r'clock_bias\(master=(\w+)\)'
    isb_pattern = r'ISB\(master=(\w+) slave=(\w+)\)'

    # Search for the pattern in the input string
    clock_match = re.search(clock_pattern, clock_string)
    isb_match = re.search(isb_pattern, isb_string)

    # If a match is found, extract the master and slave values
    if clock_match:
        master = clock_match.group(1)
        if isb_match:
            slave = isb_match.group(2)
    return master, slave
