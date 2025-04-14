class Ambiguity:
    def __init__(self, val, cov):
        self.val = val
        self.cov = cov

    def clone(self):
        obj = Ambiguity(self.val, self.cov)
        return obj




