from . import Filter


class TypeConsistencyFilter(Filter):

    def __init__(self, types):
        super().__init__()
        self.types = types

    def is_applicable(self, sat, epoch, observation):
        # return False to keep this observable
        return observation.datatype not in self.types

