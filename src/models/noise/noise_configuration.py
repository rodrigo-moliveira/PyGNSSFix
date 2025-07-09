from src.io.config import EnumNoiseProcess
from .noise_gen import *


class NoiseModel:
    # TODO: missing docstrings
    def __init__(self, name, process_str, process_noise):
        self.name = name
        self.process_enum = EnumNoiseProcess.init_model(process_str)

        if self.process_enum == EnumNoiseProcess.WHITE_NOISE:
            self.process_gen = WhiteNoiseProcess()
        elif self.process_enum == EnumNoiseProcess.RANDOM_CONSTANT:
            self.process_gen = RandomConstantProcess()
        elif self.process_enum == EnumNoiseProcess.RANDOM_WALK:
            self.process_gen = RandomWalkProcess()
        elif self.process_enum == EnumNoiseProcess.GAUSS_MARKOV:
            self.process_gen = GaussMarkovProcess()
        else:
            raise AttributeError(f"...")

        self.process_noise = process_noise

    def gen(self):
        pass

    def get_stm_entry(self):
        pass

    def get_process_noise(self):
        pass
