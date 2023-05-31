from PositioningSolver.src.gnss.data_types.ObservationData import ObservationData
from . import Functor


# Mapper class -> apply function to an input ObservationData and store the output in the output ObservationData
# This is useful for creating new ObservationData sets (for example, computing Iono Free or Smooth data)
class FunctorMapper:

    def __init__(self, functor):
        if not isinstance(functor, Functor):
            raise AttributeError(f"provided argument must be of type 'Functor'. Provided arg is {type(functor)}")
        self.functor = functor

    def apply(self, obs_data_in: ObservationData, obs_data_out: ObservationData):
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
