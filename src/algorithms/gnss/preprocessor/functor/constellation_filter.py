from src.algorithms.gnss.preprocessor.functor import Functor


class ConstellationFunctor(Functor):

    def __init__(self, constellation):
        super().__init__()
        self.constellation = constellation

    def __call__(self, obs_data_in, epoch, sat):
        if sat.sat_system == self.constellation:
            return obs_data_in.get_observables_at_epoch(epoch, sat)
        else:
            return []
