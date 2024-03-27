from src.data_mng.gnss.observation_data import ObservationData


class FilterMapper:

    def __init__(self, _filter):
        self.filter = _filter
        self.report = str()

    def apply(self, obs_data: ObservationData):
        vEpochs = list(obs_data.get_epochs())
        removed = 0
        total = 0
        sats_filter = set()

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
                for obs in v_removable:
                    obs_data.remove_observable(sat, epoch, obs.datatype)
        self.write_report(removed, total, sats_filter)

    def write_report(self, removed, total, sats_filter):
        self.report = f"The filter removed {removed/total*100:.2f}% of the data, for a total of {total} " \
                      f"observations. Satellites affected were: {sats_filter}"
