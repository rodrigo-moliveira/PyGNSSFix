""" Satellite Health Filter Module """
from . import Filter


class SatFilterHealthURA(Filter):

    def __init__(self, navigation_data, gps_ura_check, gps_ura_val, gps_health,
                 gal_sisa_check, gal_sisa_val, gal_health, log, trace_path):
        """
        Constructor for the Satellite Health and URA Filter.
        The `SatFilterHealthURA` cleans the Observation dataset by observation data of faulty satellites. The following
        checks are performed:
            * GPS URA Check: removes the GPS satellite if broadcast URA (User Range Accuracy) value is below the
                defined threshold
            * GAL SISA Check: removes the GAL satellite if broadcast SISA (Signal in Space Accuracy) value is below the
                defined threshold
            * Health Status: removes the satellite if the broadcast health status flag is faulty

        Args:
            navigation_data():
            gps_ura_check(bool): if True, the GPS URA check is performed
            gps_ura_val(float): threshold for the GPS URA
            gps_health(bool): if True, the GPS health status is evaluated
            gal_sisa_check(bool): if True, the GAL SISA check is performed
            gal_sisa_val(float): threshold for the GAL SISA
            gal_health(bool): if True, the GAL health status is evaluated
            log(logging.Logger): logger instance
            trace_path(str): path to the trace file
        """
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
        self.write_header(trace_path)

    def write_header(self, trace_path):
        if trace_path is not None:
            self.fd = open(trace_path + "/SvURAHealthFilter.txt", "w")
            self.fd.write("Epoch, Satellite, Observable, To Remove\n")

    def apply(self, sat, epoch, observation, v_removable):
        # return False to keep this observable

        flagged = False

        try:
            # get closest nav message for this satellite
            nav_message = self.navigation_data.get_closest_message(sat, epoch)
        except Exception as e:
            # no nav message available for this satellite
            if sat not in self.warned:
                self.log.warning(f"Unable to perform preprocessing filter URA/SISA and Health Check for satellite"
                                 f" {sat} due to lack of navigation data at epoch {epoch}. Exception : {e}")
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

        if self.fd is not None:
            self.fd.write(f"{epoch}, {sat}, {observation.datatype}, {flagged}\n")
