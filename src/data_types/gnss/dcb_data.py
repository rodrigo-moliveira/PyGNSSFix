class DCBData:

    def __init__(self):
        self._dcb = {"GPS": None, "GAL": None}  # DCB data from DCB files
        self._bgd = {"GPS": None, "GAL": None}  # BGD/TGD data from GAL/GPS navigation messages

    def get_dcb(self, sat, freq):
        pass

    def set_dcb(self, sat, freq, dcb):
        pass

    def get_bgd(self, sat, freq):
        pass

    def set_bgd(self, sat, freq, bgd):
        pass
