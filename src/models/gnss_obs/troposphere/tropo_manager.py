from src.errors import UnknownModel
from src.io.config.enums import EnumTropoModel, EnumTropoMask, EnumOnOff
from src.models.gnss_obs.troposphere.tropo_gpt import GPT3Tropo
from src.models.gnss_obs.troposphere.tropo_saastamoinen import SaastamoinenTropo
from src.models.gnss_obs.troposphere.mapping_function import *


class TropoManager:
    def __init__(self, config):
        tropo_model_enum = EnumTropoModel.init_model(config.get("model", "troposphere", "model"))
        tropo_mask_enum = EnumTropoMask.init_model(config.get("model", "troposphere", "mask"))
        self._estimate_tropo_wet = EnumOnOff(config.get("model", "troposphere", "estimate_tropo_wet"))

        if tropo_model_enum == EnumTropoModel.SAASTAMOINEM:
            self.tropo_model = SaastamoinenTropo()
            # when Saastamoinen is selected, the map must also be Saastamoinen
            tropo_mask_enum = EnumTropoMask.SAASTAMOINEM
        elif tropo_model_enum == EnumTropoModel.GPT3:
            from src import WORKSPACE_PATH
            from src.io.config import config_dict
            gpt_file = config_dict.get("inputs", "tropo_file")
            self.tropo_model = GPT3Tropo(WORKSPACE_PATH / f"{gpt_file}")
        elif tropo_model_enum == EnumTropoModel.DISABLED:
            self.tropo_model = None
        else:
            raise UnknownModel(f'Unknown tropo model {config.get("model", "troposphere", "model")}')

        if tropo_mask_enum == EnumTropoMask.SAASTAMOINEM:
            self.tropo_mask = SaastamoinenMap()
        elif tropo_mask_enum == EnumTropoMask.GMF:
            self.tropo_mask = GMF()
        elif tropo_mask_enum == EnumTropoMask.VMF1:
            self.tropo_mask = VMF1()
        elif tropo_mask_enum == EnumTropoMask.VMF3:
            self.tropo_mask = VMF3()
        else:
            raise UnknownModel(f'Unknown tropo mask function {config.get("model", "troposphere", "mask")}')

    def estimate_tropo(self):
        return self._estimate_tropo_wet == EnumOnOff.ENABLED

    def compute_tropo_delay(self, lat, long, height, el, epoch, state):
        if self.tropo_model is None:
            return 0.0, 0.0

        # fix height: in the estimation process, the initial state is usually [0, 0, 0] in xyz components
        # which may lead to lat,long,height=[0,0,-6378137]
        # and this produces an error in the GPT3 model, due to the bad height
        height = -height if height < -1000000 else height

        # get hydrostatic and wet delays from model
        zhd, zwd, ah, aw = self._compute_model(lat, long, height, epoch)

        if self.estimate_tropo():
            zwd = state.tropo_wet  # fetch zwd from state

        # compute mapping function coefficients
        map_hydro, map_wet = self._compute_mask(ah, aw, epoch.mjd, lat, long, height, el)

        # compute total delay
        return zhd * map_hydro + zwd * map_wet, map_wet

    def _compute_mask(self, ah, aw, mjd, lat, lon, h_ell, el):
        map_hydro = 0.0
        map_wet = 0.0
        zd = np.pi / 2.0 - el  # convert elevation to zenith angle

        if isinstance(self.tropo_mask, SaastamoinenMap):
            map_hydro, map_wet = self.tropo_mask.compute(el)
        elif isinstance(self.tropo_mask, VMF1):
            map_hydro, map_wet = self.tropo_mask.compute(ah, aw, mjd, lat, h_ell, zd)
        elif isinstance(self.tropo_mask, VMF3):
            map_hydro, map_wet = self.tropo_mask.compute(ah, aw, mjd, lat, lon, h_ell, zd)
        elif isinstance(self.tropo_mask, GMF):
            map_hydro, map_wet = self.tropo_mask.compute(mjd, lat, lon, h_ell, zd)

        return map_hydro, map_wet

    def _compute_model(self, lat, long, height, epoch):
        zhd = 0.0
        zwd = 0.0
        ah = 0.0
        aw = 0.0

        if isinstance(self.tropo_model, SaastamoinenTropo):
            doy = epoch.doy
            zhd, zwd = SaastamoinenTropo.compute(height, lat, doy)
        elif isinstance(self.tropo_model, GPT3Tropo):
            mjd = epoch.mjd
            zhd, zwd, ah, aw = self.tropo_model.compute(mjd, lat, long, height, it=0)

        return zhd, zwd, ah, aw
