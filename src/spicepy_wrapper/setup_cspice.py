import spiceypy
import os

from src.errors import CSpiceWrapperError


def create_metakernel(kernels_path, log=None):
    file_path = kernels_path / 'metakernel.tm'
    with open(file_path, 'w') as file:
        str_first = "KPL/MK \n \\begindata \n\t\t PATH_VALUES       = (\n\t\t\t'"
        files_path = str(kernels_path).replace('\\', '/') + '/'
        #files_path = str(kernels_path) + '\\'
        str_end = "'\n\t\t)\n\t\tPATH_SYMBOLS      = ('DATA')\n" \
                  "\t\tKERNELS_TO_LOAD   = (\n" \
                  "\t\t\t'$DATA\\de421.bsp',\n" \
                  "\t\t\t'$DATA\\naif0012.tls.pc',\n" \
                  "\t\t\t'$DATA\\earth_latest_high_prec.bpc',\n" \
                  "\\begintext"
        str_doc = str_first + files_path + str_end
        file.write(str_doc)
        if log is not None:
            log.info(f'CSpice Metakernel file created at {file_path}.')


def load_kernel(kernels_path, log=None):
    file_path = os.path.join(kernels_path, 'metakernel.tm')
    spiceypy.kclear()
    spiceypy.furnsh(file_path)
    if log is not None:
        log.info('Successfully loaded CSpice kernel file.')


def setup_cspice(kernels_path, log=None):
    try:
        create_metakernel(kernels_path, log)
        load_kernel(kernels_path, log)
    except Exception as e:
        raise CSpiceWrapperError(f"Error setting up CSpice: {str(e)}")
