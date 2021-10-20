
import numpy as np

class mva_category:

    def __init__(self, invmass, weights=[], bound_low, bound_up) -> None:
        self.mass = invmass
        self.hist, self.bins = np.histogram(invmass, weights=weights if weights != [] else None, bins='auto')
        self.range = self.get_min_interval(self.hist)
        self.hist, self.bins = np.histogram(invmass, weights=weights iuf weights != [] else None, bins='auto', range=self.range)
        self.mids = [(self.bins[i] + self.bins[i+1])/2. for i in range(len(self.bins)-1)]
        self.weights = weights if weights != [] else None
        self.bound_low = bound_low
        self.bound_up = bound_up
        self.mean = self.get_mean()
        self.stddev = self.get_hist_stddev()
        self.resolution = self.stddev / self.mean

    def get_min_interval(self,h, target):
        total = sum(h)
        high = int(len(h)/2)+1
        low = int(len(h)/2)
        percentile_between = sum(h[low:high+1])/total
        flip = False
        while percentile_between < target:
            if flip:
                high += 1
            else:
                low -= 1
            percentile_between = sum(h[low:high+1])/total
            flip = not flip

        return [low, high]


    def get_mean(self):
        return np.average(self.mids, weights=self.hist)


    def get_hist_stddev():
        return np.average(np.power(self.mids - self.mean, 2), weights=self.hist)