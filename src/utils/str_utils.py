import re

def replace_whitespace_with_underscore(s):
    # Replace one or more whitespace characters with a single underscore
    return re.sub(r'\s+', '_', s)
