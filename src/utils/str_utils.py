import re


def replace_whitespace_with_underscore(s):
    # Replace one or more whitespace characters with a single underscore
    return re.sub(r'\s+', '_', s)


def extract_constellations(clock_string, isb_string):
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

