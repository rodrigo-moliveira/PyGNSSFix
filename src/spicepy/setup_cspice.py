import spiceypy
import os


def create_metakernel(kernels_path):
    file_path = os.path.join(kernels_path, 'metakernel.tm')
    with open(file_path, 'w') as file:
        str_first = "KPL/MK \n \\begindata \n\t\t PATH_VALUES       = (\n\t\t\t'"
        files_path = kernels_path.replace('\\', '/') + '/'
        str_end = "'\n\t\t)\n\t\tPATH_SYMBOLS      = ('DATA')\n\t\tKERNELS_TO_LOAD   = (\n\t\t\t'$DATA/de430.bsp'," \
                  "\n\t\t\t'$DATA/naif0012.tls.pc',\n\\begintext"
        str_doc = str_first + files_path + str_end
        file.write(str_doc)


def load_kernel(kernels_path):
    file_path = os.path.join(kernels_path, 'metakernel.tm')
    spiceypy.kclear()
    print(f'Loading CSpice kernel file {file_path}...')
    spiceypy.furnsh(file_path)
    print('CSpice library successfully installed.')


def setup_cspice():
    _kernels_path = 'C:\\Users\\rooo\\Documents\\other_projects\\PyGNSSFix\\workspace\\cspice'
    create_metakernel(_kernels_path)
    load_kernel(_kernels_path)