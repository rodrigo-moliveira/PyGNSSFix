""" This module contains the functions to set up the spicepy library
(pythonic wrapper of the CSpice library). """

import spiceypy
import os

from src.errors import CSpiceWrapperError


def create_metakernel(kernels_path, kernel_list, log=None):
    """ Create the CSpice metakernel file (`metakernel.tm`).

    Args:
        kernels_path (str): Path to the kernels directory.
        kernel_list (list): List of kernel files to load.
        log (logging.Logger): Logger object to log messages.
    """
    file_path = kernels_path + '\\metakernel.tm'
    with open(file_path, 'w') as file:
        str_first = "KPL/MK \n \\begindata \n\t\t PATH_VALUES       = (\n\t\t\t'"
        files_path = str(kernels_path).replace('\\', '/') + '/'
        str_end = "'\n\t\t)\n\t\tPATH_SYMBOLS      = ('DATA')\n" \
                  "\t\tKERNELS_TO_LOAD   = (\n"
        for kernel in kernel_list:
            str_end += f"\t\t\t'$DATA\\{kernel}',\n"
        str_end += "\t\t)\n\\begintext"

        str_doc = str_first + files_path + str_end
        file.write(str_doc)

    # Log if provided
    if log is not None:
        log.info(f'CSpice Metakernel file created at {file_path}.')
    else:
        print(f'CSpice Metakernel file created at {file_path}.')


def load_kernel(kernels_path, log=None):
    """ Load the CSpice kernel file.
    Calls the `furnsh` function from the `spiceypy` library to load the kernel file.

    Args:
        kernels_path (str): Path to the kernels directory.
        log (logging.Logger): Logger object to log messages.
    """
    file_path = os.path.join(kernels_path, 'metakernel.tm')
    spiceypy.kclear()
    spiceypy.furnsh(file_path)
    if log is not None:
        log.info('Successfully loaded CSpice kernel file.')
    else:
        print('Successfully loaded CSpice kernel file.')


def setup_cspice(kernels_path, kernel_list, log=None):
    """ Set up the `spiceypy` library (wrapper of CSpice).
     Creates and loads the required kernel files.

    Args:
        kernels_path (str): Path to the kernel files directory.
        kernel_list (list): List of kernel files to load.
        log (logging.Logger): Logger object to log messages.
    """
    try:
        # get relative path
        current_dir = os.getcwd()
        relative_path = os.path.relpath(kernels_path, current_dir)

        # create metakernel file
        create_metakernel(relative_path, kernel_list, log)

        # load kernel
        load_kernel(relative_path, log)
    except Exception as e:
        raise CSpiceWrapperError(f"Error setting up CSpice: {str(e)}")
