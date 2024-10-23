""" Filter Mapper (Transformer) Module """
from src.data_mng.gnss.observation_data import ObservationData
from .filter import Filter


class FilterMapper:
    """
    Manages the execution of the filter algorithms and effectively removes the flagged data from the
    ObservationData instance.

    Attributes:
        filter(Filter): Filter to be applied
        report(str): after the Filter is applied, a small report is written for logging purposes.
    """

    def __init__(self, filter_instance: Filter):
        """
        Constructor of the FilterMapper.

        Args:
            filter_instance(Filter): instance of the `Filter` class
        """
        self.filter = filter_instance
        self.report = str()

    def apply(self, obs_data: ObservationData):
        """
        Applies the filter to the provided `ObservationData` instance.
        Iterates over all epochs, satellites and observables and calls the `filter.is_applicable` and `filter.apply`
        methods (must be implemented by the inherited Filter class)

        Args:
            obs_data(ObservationData): input observation dataset to be filtered (data is effectively removed from the
                input dataset)
        """
        vEpochs = list(obs_data.get_epochs())
        removed = 0
        total = 0
        sats_filter = set()
        removed_epochs = set()
        total_epochs = len(vEpochs)

        # epoch loop
        for epoch in vEpochs:
            epoch_data = obs_data.get_epoch_data(epoch)

            # get available satellites
            vSats = list(epoch_data.get_satellites())

            # satellite loop
            for sat in vSats:

                # list of observables to remove
                v_removable = []
                v_observables = obs_data.get_observables_at_epoch(epoch, sat)
                total += len(v_observables)

                for obs in v_observables:
                    if self.filter.is_applicable(sat, epoch, obs, obs_list=v_observables):
                        self.filter.apply(sat, epoch, obs, v_removable)

                # remove observables
                removed += len(v_removable)
                if len(v_removable) > 0:
                    sats_filter.add(sat)
                    removed_epochs.add(epoch)
                for obs in v_removable:
                    obs_data.remove_observable(sat, epoch, obs.datatype)
        self.write_report(removed, total, sats_filter, total_epochs, len(removed_epochs))

    def write_report(self, removed, total, sats_filter, total_epochs, removed_epochs):
        self.filter.close_file()
        self.report = f"The filter removed {removed/total*100:.2f}% of the data, for a total of {total} " \
                      f"observations. Satellites affected were: {sats_filter}. Percentage of epochs affected is " \
                      f"{removed_epochs/total_epochs*100:.2f}%."
