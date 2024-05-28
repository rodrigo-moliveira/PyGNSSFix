""" Functor Mapper (Transformer) Module """

from src.data_mng.gnss.observation_data import ObservationData
from . import Functor


class FunctorMapper:
    """
    Manages the execution of the functor algorithms, that apply a given transformation function from the input
    `ObservationData` instance to the output `ObservationData` instance.

    For example, for the computation of the Smooth Observation Data
    """

    def __init__(self, functor):
        """
        Constructor of the FunctorMapper

        Parameters:
            functor(Functor): instance of the `Functor` class
        """
        if not isinstance(functor, Functor):
            raise AttributeError(f"provided argument must be of type 'Functor'. Provided arg is {type(functor)}")
        self.functor = functor

    def apply(self, obs_data_in: ObservationData, obs_data_out: ObservationData):
        """
        Applies the functor to the input Observation data and stores the result in the output Observation data
        Iterates over all epochs and satellites and calls the `functor.__call__` method of the functor
        (must be implemented by the inherited Functor class)

        Parameters:
            obs_data_in(ObservationData): input observation dataset (read-only)
            obs_data_out(ObservationData): output observation dataset
        """
        vEpochs = obs_data_in.get_epochs()

        # epoch loop
        for epoch in vEpochs:
            epoch_data = obs_data_in.get_epoch_data(epoch)

            # get available satellites
            vSats = epoch_data.get_satellites()

            # satellite loop
            for sat in vSats:

                vObs = self.functor(obs_data_in, epoch, sat)

                # set observables
                for obs in vObs:
                    obs_data_out.set_observation(epoch, sat, obs)
