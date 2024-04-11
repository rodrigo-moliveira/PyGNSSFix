from src.io.config import EnumIonoModel, config_dict as config
from src.errors import ConfigError
from .iono_klobuchar import IonoKlobuchar
from .iono_ntcmg import NTCMG
from src.common_log import get_logger, MODEL_LOG


# NOTE: a cache interface may be implemented, for example to re-use the iono computation for the previous
# epoch, for example (if the receiver position does not change much)

class IonoManager:
    """
    Class with Management of a-priori Ionosphere models
    Implemented models:
        * Klobuchar (recommended for GPS satellites)
        * NTCM-G (recommended for GAL satellites)

    Note that each constellation provider developed its own iono model, however the software allows
    to use any model desired (for example, we can run GPS with NTCM-G or GAL with Klobuchar without
    any problem).
    """
    def __init__(self, constellation):
        """
        Constructor of the Iono Manager
        Args:
            constellation(str): input constellation to associate to this IonoManager object
        """
        iono_model_enum = EnumIonoModel.init_model(config.get("model", constellation, "ionosphere"))

        if iono_model_enum == EnumIonoModel.KLOBUCHAR:
            self.iono_model = IonoKlobuchar()
        elif iono_model_enum == EnumIonoModel.NTCMG:
            self.iono_model = NTCMG()
        else:
            raise ConfigError(f'Unknown iono model {iono_model_enum}')

    def compute_iono_delay(self, epoch, iono_corrections, sat, user_lat, user_long, sv_el, sv_az, freq):
        """
        Compute iono delay in meters for the corresponding a-priori model selected
        Args:
            epoch(src.data_types.date.date.Epoch): required epoch to compute the iono delay
            iono_corrections(dict): dictionary with iono corrections read from RINEX NAV header
            sat(str): satellite
            user_lat(float): user latitude in [rad]
            user_long(float): user longitude in [rad]
            sv_el(float): satellite elevation in [rad]
            sv_az(float): satellite azimuth in [rad]
            freq(src.data_types.gnss.data_type.DataType): Frequency band required for the iono delay
        Returns:
            float: computed iono delay in meters for the provided frequency band
        """
        iono = 0.0
        if isinstance(self.iono_model, IonoKlobuchar):
            epoch_gps = epoch.change_scale("GPST")
            gps_sow = epoch_gps.gnss_time[1]
            try:
                alfa = iono_corrections["GPSA"]
                beta = iono_corrections["GPSB"]
                if len(alfa) != 4 or len(beta) != 4:
                    raise AttributeError(f"Alfa and Beta iono vectors must have dimension 4. "
                                         f"Provided vectors: Alfa={alfa}. Beta={beta}")
                iono = self.iono_model.compute(gps_sow, alfa, beta, user_lat, user_long, sv_el, sv_az, freq)
            except Exception as e:
                log = get_logger(MODEL_LOG)
                log.warning(f"Problem at computing iono a-priori model Klobuchar for {sat} at {epoch}: {e}")

        elif isinstance(self.iono_model, NTCMG):
            ut1 = epoch.change_scale("UT1").datetime
            try:
                gal_param = iono_corrections["GAL"]
                if len(gal_param) != 3:
                    raise AttributeError(f"Galileo iono corrections vector must have dimension 3. "
                                         f"Provided vector: GAL iono={gal_param}")
                iono = self.iono_model.compute(ut1, gal_param, user_lat, user_long, sv_el, sv_az, freq)
            except Exception as e:
                log = get_logger(MODEL_LOG)
                log.warning(f"Problem at computing iono a-priori model NTCM-G for {sat} at {epoch}: {e}")
        return iono
