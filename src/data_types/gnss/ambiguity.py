class Ambiguity:
    def __init__(self, val, cov):
        self.val = val
        self.cov = cov

    def clone(self):
        obj = Ambiguity(self.val, self.cov)
        return obj

    def __str__(self):
        return f"Ambiguity(val={self.val}, cov={self.cov})"

    def __repr__(self):
        return str(self)




