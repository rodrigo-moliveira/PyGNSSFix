
_cache = {}
class ITRF_Transformation:
    def __init__(self, in_frame, out_frame, helmert_data):
        self.in_frame = in_frame
        self.out_frame = out_frame
        self.helmert_data = helmert_data

