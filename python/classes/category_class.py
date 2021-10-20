

class mva_category:

    def __init__(self, invmass, weights=[], bound_low, bound_up) -> None:
        self.mass = invmass
        self.weights = weights,
        self.bound_low = bound_low
        self.bound_up = bound_up
        self.resolution