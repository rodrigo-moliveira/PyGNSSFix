""" GNSS Algorithm Manager Module """

import traceback
import numpy as np

from src.io.config import config_dict, EnumAlgorithmINS
from src.errors import ConfigError
from src.common_log import MAIN_LOG, get_logger
from src.data_mng.ins.ins_data_mng import InsDataManager

#from .solver.gnss_solver import GnssSolver
#from .preprocessor import PreprocessorManager
#from ...data_mng.gnss.gnss_data_mng import GnssDataManager


class InsAlgorithmManager:
    """
    TODO: to be updated
    GNSS Algorithm Manager.

    This class executes the GNSS PNT solver algorithm. It is structured in the following tasks:
        * Constructs the :py:class:`GnssDataManager` instance
        * Reads all the input data (:py:meth:`GnssDataManager.read_inputs`)
        * Launches the main PNT algorithm (function `compute`), which performs the following sub-tasks:
            * launches the :py:class:`PreprocessorManager` algorithm
            * launches the :py:class:`GnssSolver` algorithm
            * compute DOPs (see :py:meth:`compute_dop`)

        Attributes:
            data_dir(str): absolute path where the output data is stored
            data_manager(GnssDataManager): instance of the GNSS Data Manager
            main_log(logging.Logger): logger instance
    """

    def __init__(self, data_dir):
        # create data members
        self.data_manager = InsDataManager()
        self.main_log = None
        self.data_dir = data_dir

    def run(self):
        """ Main function that executes the INS Algorithm """
        self.main_log = get_logger(MAIN_LOG)
        self.main_log.info("Starting INS Algorithm Manager")
        ins_alg = config_dict.get('ins_alg')

        if ins_alg not in (EnumAlgorithmINS.SENSOR_EMULATOR, EnumAlgorithmINS.INS_INTEGRATION, EnumAlgorithmINS.LOOSELY_COUPLED):
            raise ConfigError(f"Selected Model {ins_alg} not valid. Available options are "
                              f"Sensor Emulator, INS Integration, Loosely Coupled")
        self.main_log.info(f"Running INS algorithm {ins_alg}")

        # Input Reader Module
        try:
            self.main_log.info(f"Starting Input Reader Module...")
            self.data_manager.read_inputs(ins_alg, f"{self.data_dir}\\trace")
        except Exception as e:
            self.main_log.error(f"Stopping execution of program due to error in execution of Input Reader Module: {e}")
            print(traceback.format_exc())
            exit(-1)

        # Main Algorithm Module
        #try:
        #    self.main_log.info(f"Starting Main Algorithm Module...")
        #    self._compute(self.data_manager, f"{self.data_dir}\\trace")
        #except Exception as e:
        #    self.main_log.error(f"Stopping execution of program due to error in execution of Main Algorithm "
        #                        f"Module: {e}")
        #    print(traceback.format_exc())
        #    exit(-1)

        # Output Writer Module
        #try:
        #    self.main_log.info(f"Starting Output Writer Module...")
        #    self.data_manager.save_data(f"{self.data_dir}\\output")
        #except Exception as e:
        #    self.main_log.error(f"Stopping execution of program due to error in execution of Output Writer Module: "
        #                        f"{str(e)}")
        #    print(traceback.format_exc())
        #    exit(-1)

        self.main_log.info(f"Successfully executed INS algorithm {ins_alg}")


    #def _compute(self, data_manager, trace_path):
    #    """ Internal function for the computation of the PNT solver task """
    #    # perform pre-processing here
    #    self.main_log.info(f"Starting Preprocessor Module")
    #    preprocessor = PreprocessorManager(trace_path, data_manager)
    #    preprocessor.compute()

        # run estimation algorithm
    #    self.main_log.info(f"Running estimation algorithm...")
    #    solver = GnssSolver(data_manager, trace_path)
    #    solver.solve()

    #    data_manager.add_data("nav_solution", solver.solution)

