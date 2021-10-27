
import numpy as np
from scipy.stats import iqr

class mva_category:

    def __init__(self, invmass, weights, bound_low, bound_up) -> None:
        if bound_low == bound_up:
            raise Exception("no events, bounds are equal")
        self.mass = np.array(invmass)
        self.weights = np.array(weights)
        self.range = self.weighted_quantile(0.683)
        mask = np.logical_and(self.range[0] <= invmass, invmass <= self.range[1])
        self.mass = invmass[mask]
        self.weights = self.weights[mask]
        self.lower = bound_low
        self.upper = bound_up
        self.mean = self.get_mean()
        self.stddev = self.get_weighted_stddev()
        self.resolution = self.stddev / self.mean

    def weighted_quantile(self, quantiles):
        """ calculates unbinned interval containing target percent of events """
        quantiles = [0.5-(quantiles/2), 0.5+(quantiles/2)]
        sorter = np.argsort(self.mass)
        values = self.mass[sorter]
        sample_weight = self.weights[sorter]
        weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
        weighted_quantiles /= np.sum(sample_weight)
        return np.interp(quantiles, weighted_quantiles, values)

    def get_mean(self):
        return np.average(self.mass, weights=self.weights)

    def get_weighted_stddev(self):
        return np.sqrt(np.average(np.power((self.mass - self.mean),2), weights=self.weights))
