""" Module with the base Functor Class """


class Functor:
    """
    Base class for Functor classes.
    A functor is an algorithm that transforms an `ObservationData` object to another, implementing a function
    with the observables.

    Examples of functors are:
        * IonoFree functor
        * Smooth functor

    Functor classes inherited from this must implement the ``call`` method.
    """

    def __init__(self):
        pass

    def __call__(self, obs_data_in, epoch, sat):
        """
        Implementation of the functor algorithm - `call` method.

        Args:
            obs_data_in(src.data_mng.gnss.observation_data.ObservationData): `ObservationData` instance
            epoch(src.data_types.date.Epoch): epoch under evaluation
            sat(src.data_types.gnss.Satellite): satellite under evaluation

        Returns:
            list: list with processed observables to insert in the new `ObservationData` instance
        """
        raise NotImplemented(f"Functor not yet implement. It must define a `call` method")
