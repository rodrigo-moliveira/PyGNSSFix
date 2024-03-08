from . import Filter


class SatFilterHealthURA(Filter):

    def __init__(self, navigation_data, gps_ura_check, gps_ura_val, gps_health,
                 gal_sisa_check, gal_sisa_val, gal_health, log):
        super().__init__()
        self.navigation_data = navigation_data
        self.gps_ura_check = gps_ura_check
        self.gps_ura_val = gps_ura_val
        self.gps_health = gps_health
        self.gal_sisa_check = gal_sisa_check
        self.gal_sisa_val = gal_sisa_val
        self.gal_health = gal_health
        self.log = log
        self.warned = []

    def apply(self, sat, epoch, observation, v_removable):
        # return False to keep this observable

        flagged = False

        try:
            # get closest nav message for this satellite
            nav_message = self.navigation_data.get_closest_message(sat, epoch)
        except:
            # no nav message available for this satellite
            if sat not in self.warned:
                self.log.warn(f"Unable to perform preprocessing filter URA/SISA and Health Check for satellite"
                              f" {sat} due to lack of navigation data at epoch {epoch}")
                self.warned.append(sat)
            return False

        # GPS checks
        if sat.sat_system == "GPS":
            if self.gps_ura_check and nav_message.SV_URA > self.gps_ura_val:
                self.log.debug(f"Satellite {sat} is being discarded at epoch {str(epoch)} due to high "
                               f"URA value ({nav_message.SV_URA})m compared to threshold {self.gps_ura_val}m")
                flagged = True

            # Health check
            if self.gps_health and nav_message.SV_health != 0:
                self.log.debug(f"Satellite {sat} is being discarded at epoch {str(epoch)} due to bad "
                               f"health flag ({nav_message.SV_health}) in the navigation message")
                flagged = True

        # GAL checks
        elif sat.sat_system == "GAL":
            if self.gal_sisa_check and nav_message.SISA > self.gal_sisa_val:
                self.log.debug(f"Satellite {sat} is being discarded at epoch {str(epoch)} due to high "
                               f"SISA value ({nav_message.SISA})m compared to threshold {self.gal_sisa_val}m")
                flagged = True

        if flagged:
            v_removable.append(observation)
