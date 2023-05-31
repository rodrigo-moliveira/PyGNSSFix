from .filter import Filter


class RateDowngradeFilter(Filter):

    def __init__(self, rate_out: float, first_epoch):
        super().__init__()
        self.rate_out = rate_out
        self.first_epoch = first_epoch

    def is_applicable(self, sat, epoch, observation):
        # return True to remove observables, return False to keep observable
        return (epoch - self.first_epoch) % self.rate_out != 0
