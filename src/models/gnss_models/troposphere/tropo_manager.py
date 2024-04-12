import numpy as np

from src.errors import ConfigError
from src.io.config.enums import EnumTropoModel, EnumTropoMap, EnumOnOff
from src.models.gnss_models.troposphere.gpt import *
from src.models.gnss_models.troposphere.tropo_saastamoinen import SaastamoinenTropo
from src.io.config.config import config_dict as config


class TropoManager:
    """
    Class with Management of a-priori Troposphere models

    Implemented Models:
        * Saastamoinen
        * GPT3

    The user can also choose the map function. The following options are available:
        * Saastamoinen
        * GMF
        * VMF1
        * VMF3

    The user may also choose if the wet component of the troposphere can be estimated in the algorithm
    or the a-priori one is used instead.
    """
    def __init__(self):
        # get user configurations
        tropo_model_enum = EnumTropoModel.init_model(config.get("model", "troposphere", "model"))
        tropo_map_enum = EnumTropoMap.init_model(config.get("model", "troposphere", "map"))
        self._estimate_tropo_wet = EnumOnOff(config.get("model", "troposphere", "estimate_tropo_wet"))

        # verifications
        if tropo_model_enum == EnumTropoModel.SAASTAMOINEM:
            self.tropo_model = SaastamoinenTropo()
            # when Saastamoinen is selected, the map must also be Saastamoinen
            tropo_map_enum = EnumTropoMap.SAASTAMOINEM
        elif tropo_model_enum == EnumTropoModel.GPT3:
            from src import WORKSPACE_PATH
            from src.io.config import config_dict
            gpt_file = config_dict.get("inputs", "tropo_file")
            self.tropo_model = GPT3Tropo(WORKSPACE_PATH / f"{gpt_file}")
        elif tropo_model_enum == EnumTropoModel.DISABLED:
            self.tropo_model = None
        else:
            raise ConfigError(f'Unknown tropo model {config.get("model", "troposphere", "model")}')

        if tropo_map_enum == EnumTropoMap.SAASTAMOINEM:
            self.tropo_map = SaastamoinenMap()
        elif tropo_map_enum == EnumTropoMap.GMF:
            self.tropo_map = GMF()
        elif tropo_map_enum == EnumTropoMap.VMF1:
            self.tropo_map = VMF1()
        elif tropo_map_enum == EnumTropoMap.VMF3:
            self.tropo_map = VMF3()
        else:
            raise ConfigError(f'Unknown tropo map function {config.get("model", "troposphere", "map")}')

    def estimate_tropo(self):
        """
        Returns:
            bool: True if wet component of tropo is estimated in the filter algorithm
        """
        return self._estimate_tropo_wet == EnumOnOff.ENABLED

    def compute_tropo_delay(self, lat, long, height, el, epoch, state_zwd):
        """
        Compute tropospheric delay.
        If the estimation of the zenith wet delay (zwd) component is active, then the `state_zwd` parameter
        (estimated by the filter) is used. Otherwise, the zwd from the model is used instead

        Args:
            lat(float): user latitude in [rad]
            long(float): user longitude in [rad]
            height(float): user height in [m]
            el(float): satellite elevation from user position in [rad]
            epoch(src.data_types.date.date.Epoch): epoch to compute the tropo delay
            state_zwd(float or None): state estimated zenith wet delay (wzd) or None

        Returns:
            tuple: A tuple with tropo_delay in meters and wet map function is returned as (tropo_delay, map_wet)
        """
        if self.tropo_model is None:
            return 0.0, 0.0

        # fix height: in the estimation process, the initial state is usually [0, 0, 0] in xyz components
        # which may lead to lat,long,height=[0,0,-6378137]
        # and this produces an error in the GPT3 model, due to the bad height
        height = -height if height < -1000000 else height

        # get hydrostatic and wet delays from model
        zhd, zwd, ah, aw = self._compute_model(lat, long, height, epoch)

        if self.estimate_tropo():
            # when tropo estimation is on, the model zwd is ignored, and the one from the state is used instead
            zwd = state_zwd

        # compute mapping function coefficients
        map_hydro, map_wet = self._compute_map(ah, aw, epoch.mjd, lat, long, height, el)

        # compute total delay
        return zhd * map_hydro + zwd * map_wet, map_wet

    def _compute_map(self, ah, aw, mjd, lat, lon, h_ell, el):
        """Compute hydrostatic and wet mapping functions
        """
        map_hydro = 0.0
        map_wet = 0.0
        zd = np.pi / 2.0 - el  # convert elevation to zenith angle

        if isinstance(self.tropo_map, SaastamoinenMap):
            map_hydro, map_wet = self.tropo_map.compute(el)
        elif isinstance(self.tropo_map, VMF1):
            map_hydro, map_wet = self.tropo_map.compute(ah, aw, mjd, lat, h_ell, zd)
        elif isinstance(self.tropo_map, VMF3):
            map_hydro, map_wet = self.tropo_map.compute(ah, aw, mjd, lat, lon, h_ell, zd)
        elif isinstance(self.tropo_map, GMF):
            map_hydro, map_wet = self.tropo_map.compute(mjd, lat, lon, h_ell, zd)

        return map_hydro, map_wet

    def _compute_model(self, lat, long, height, epoch):
        """Compute hydrostatic and wet zenith delays from the a-priori models
        """
        zhd = 0.0
        zwd = 0.0
        ah = 0.0
        aw = 0.0

        if isinstance(self.tropo_model, SaastamoinenTropo):
            doy = epoch.doy
            zhd, zwd = self.tropo_model.compute(height, lat, doy)
        elif isinstance(self.tropo_model, GPT3Tropo):
            mjd = epoch.mjd
            zhd, zwd, ah, aw = self.tropo_model.compute(mjd, lat, long, height, it=0)
        return zhd, zwd, ah, aw
