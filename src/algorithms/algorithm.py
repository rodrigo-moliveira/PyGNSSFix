class Algorithm:
    def __init__(self):
        self.inputs = []
        self.outputs = []
        self.results = []

    def __str__(self):
        return "Algorithm(Unknown Algorithm)"

    def compute(self, *args):
        pass

    def get_results(self, *args):
        return self.results
