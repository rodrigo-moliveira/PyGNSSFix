from PositioningSolver.src.gnss.data_types.ObservationData import ObservationData


class FilterMapper:

    def __init__(self, _filter):
        self.filter = _filter

    def apply(self, obs_data: ObservationData):
        vEpochs = list(obs_data.get_epochs())

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

                for obs in v_observables:
                    if self.filter.is_applicable(sat, epoch, obs):
                        self.filter.apply(sat, epoch, obs, v_removable)

                # remove observables
                for obs in v_removable:
                    obs_data.remove_observable(sat, epoch, obs.datatype)
