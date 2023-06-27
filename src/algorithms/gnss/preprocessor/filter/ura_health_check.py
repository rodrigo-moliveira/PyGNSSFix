from . import Filter


class SatFilterHealthURA(Filter):

    def __init__(self, navigation_data, ura_check, ura_threshold, health_check, log):
        super().__init__()
        self.navigation_data = navigation_data
        self.ura_check = ura_check
        self.ura_threshold = ura_threshold
        self.health_check = health_check
        self.log = log
        self.warned = []

    def apply(self, sat, epoch, observation, v_removable):
        # return False to keep this observable

        flagged = False

        try:
            # get closest nav message for this satellite
            nav_message = self.navigation_data.get_sat_data_for_epoch(sat, epoch)
        except:
            # no nav message available for this satellite
            if sat not in self.warned:
                self.log.warn(f"Unable to perform preprocessing filter SV URA and Health Check for satellite"
                              f" {sat} due to lack of navigation data")
                self.warned.append(sat)
            return False

        # URA check
        if self.ura_check and nav_message.SV_URA > self.ura_threshold:
            self.log.debug(f"Satellite {sat} is being discarded at epoch {str(epoch)} due to high "
                           f"URA value ({nav_message.SV_URA}) compared to threshold {self.ura_threshold}")
            flagged = True

        # Health check
        if self.health_check and nav_message.SV_health != 0:
            self.log.debug(f"Satellite {sat} is being discarded at epoch {str(epoch)} due to bad "
                           f"health flag ({nav_message.SV_health}) in the navigation message")
            flagged = True

        if flagged:
            v_removable.append(observation)
