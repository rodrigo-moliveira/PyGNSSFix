from src.errors import UnknownModel
from src.io.config.enums import EnumTropoModel, EnumTropoMask, EnumOnOff
from src.models.gnss_obs.troposphere.tropo_gpt import GPT3Tropo
from src.models.gnss_obs.troposphere.tropo_saastamoinen import SaastamoinenTropo
from src.models.gnss_obs.troposphere.mapping_function import *


class TropoManager:
    def __init__(self, config):
        tropo_model_enum = EnumTropoModel.init_model(config.get("model", "troposphere", "model"))
        tropo_mask_enum = EnumTropoMask.init_model(config.get("model", "troposphere", "mask"))
        self.estimate_tropo_wet = EnumOnOff(config.get("model", "troposphere", "estimate_tropo_wet"))

        if tropo_model_enum == EnumTropoModel.SAASTAMOINEM:
            self.tropo_model = SaastamoinenTropo()
        elif tropo_model_enum == EnumTropoModel.GPT3:
            self.tropo_model = GPT3Tropo()
        elif tropo_model_enum == EnumTropoModel.DISABLED:
            self.tropo_model = None
        else:
            raise UnknownModel("")

        if tropo_mask_enum == EnumTropoMask.SAASTAMOINEM:
            self.tropo_mask = SaastamoinenMap()
        elif tropo_mask_enum == EnumTropoMask.GMF:
            self.tropo_mask = GMF()
        elif tropo_mask_enum == EnumTropoMask.VMF1:
            self.tropo_mask = VMF1()
        elif tropo_mask_enum == EnumTropoMask.VMF3:
            self.tropo_mask = VMF3()
        else:
            raise UnknownModel("")

    def compute_tropo_delay(self):
        # get hydrostatic and wet delays from model
        zhd, zwd = self.compute_model()

        # compute mapping function coefficients
        map_hydro, map_wet = self.compute_mask()

        # compute total delay
        return zhd * map_hydro + zwd * map_wet

    def compute_mask(self, ah, aw, mjd, lat, lon, h_ell, zd):
        map_hydro = 0.0
        map_wet = 0.0

        if isinstance(self.tropo_mask, SaastamoinenMap):
            pass
        elif isinstance(self.tropo_mask, VMF1):
            map_hydro, map_wet = self.tropo_mask.compute(ah, aw, mjd, lat, h_ell, zd)
        elif isinstance(self.tropo_mask, VMF3):
            map_hydro, map_wet = self.tropo_mask.compute(ah, aw, mjd, lat, lon, h_ell, zd)
        elif isinstance(self.tropo_mask, GMF):
            map_hydro, map_wet = self.tropo_mask.compute(mjd, lat, lon, h_ell, zd)

        return map_hydro, map_wet

    def compute_model(self):
        zhd = 0.0
        zwd = 0.0

        if self.tropo_model is None:
            pass
        elif isinstance(self.tropo_model, SaastamoinenTropo):
            zhd, zwd = SaastamoinenTropo.compute(h, lat, doy)
        elif isinstance(self.tropo_model, GPT3Tropo):
            zhd, zwd = self.tropo_model.compute(mjd, lat, lon, h_ell, it=0)

        return zhd, zwd
    