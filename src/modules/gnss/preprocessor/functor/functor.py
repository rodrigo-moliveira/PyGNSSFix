# base functor class

# Functor classes inherited from this must implement the ``call`` method
class Functor:
    def __init__(self):
        pass

    def __call__(self, obs_data_in, epoch, sat):
        raise NotImplemented(f"Functor not yet implement. It must define a ``call`` method")
