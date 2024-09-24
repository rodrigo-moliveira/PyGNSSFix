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
