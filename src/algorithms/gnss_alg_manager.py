import os
import time

from src import PROJECT_PATH
from src.data_mng.gnss_data_mng import GnssDataManager


class GnssAlgorithmManager:

    def __init__(self, algorithm):
        self.data_manager = GnssDataManager()
        self.algorithm = algorithm

    def read_inputs(self):

        # read navigation data
        nav_data = None
        self.data_manager.add_data("navigation_data", nav_data)

        # read observation data
        obs_data = None
        self.data_manager.add_data("observation_data", obs_data)

        # ... add more here

    def run(self):

        # fetch input variables to this algorithm
        input_names = self.algorithm.inputs
        inputs = []

        for _in_name in input_names:
            _in = self.data_manager.get_data(_in_name)
            inputs.append(_in)

        self.algorithm.compute(*inputs)

        # get results and add them to the data manager
        _results = self.algorithm.get_results()
        for _data, _name in zip(_results, self.algorithm.outputs):
            if _data is not None:
                self.data_manager.add_data(_name, _data)

    def results(self, data_dir=None, performance=False, plot=False, separate_axis=False):
        #### check data dir

        if data_dir is not None:  # data_dir specified, meaning to save .csv files
            data_dir = self._check_data_dir(data_dir)

            # save data files
            self.data_manager.save_data(data_dir)

        # GnssQualityManager.process(self.data_manager, data_dir, self.algorithm.name, performance, plot, separate_axis)

    def _check_data_dir(self, data_dir):
        """
        check if data_dir is a valid dir. If not, use the default dir.
        check if the data_dir exists. If not, create it.
        Args:
            data_dir: all generated files are saved in data_dir
        Returns:
            data_dir: valid data dir.
        """
        # check data dir
        # data_dir is not specified, automatically create one
        if data_dir == '':
            data_dir = PROJECT_PATH
            if data_dir[-1] != '//':
                data_dir = data_dir + '//'
            data_dir = data_dir + time.strftime('%Y-%m-%d-%H-%M-%S', time.localtime()) + '//'
            data_dir = os.path.abspath(data_dir)
            print("creating dir ", data_dir)
        # create data dir
        if not os.path.exists(data_dir):
            try:
                data_dir = os.path.abspath(data_dir)
                os.makedirs(data_dir)
            except:
                raise IOError(f"Cannot create dir: {data_dir}")
        return data_dir
