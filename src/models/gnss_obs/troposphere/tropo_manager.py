from src.io.config.enums import EnumTropoModel, EnumTropoMask, EnumOnOff
from src.models.gnss_obs.troposphere.tropo_gpt import GPT3
from src.models.gnss_obs.troposphere.tropo_saastamoinen import Saastamoinen


class TropoManager:
    def __init__(self, config):
        tropo_model_enum = EnumTropoModel.init_model(config.get("model", "troposphere", "model"))
        self.tropo_mask = EnumTropoMask.init_model(config.get("model", "troposphere", "mask"))
        self.estimate_tropo_wet = EnumOnOff(config.get("model", "troposphere", "estimate_tropo_wet"))

        if tropo_model_enum == EnumTropoModel.SAASTAMOINEM:
            self.tropo_model = Saastamoinen()
        elif tropo_model_enum == EnumTropoModel.GPT3:
            self.tropo_model = GPT3()
        elif tropo_model_enum == EnumTropoModel.DISABLED:
            self.tropo_model = None

        if self.tropo_mask == EnumTropoMask.SAASTAMOINEM:
            pass
        elif self.tropo_mask == EnumTropoMask.GMF:
            pass
        elif self.tropo_mask == EnumTropoMask.VMF1:
            pass
        elif self.tropo_mask == EnumTropoMask.VMF3:
            pass

    def compute_tropo_delay(self):
        pass

    def compute_mask(self):
        pass

    def compute_model(self):
        pass
    