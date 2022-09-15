
from matplotlib.pyplot import xcorr
import numpy as np
from scipy.stats import moment


class mva_category:
    """ class to describe diphoton mva categories """

    def __init__(self, invmass, weights, is_signal) -> None:

        # evaluate range of 68.3% min interval
        self.range = self.weighted_quantile(invmass, weights, 0.683)
        if len(self.range) == 0 or sum(weights) == 0:
            self.s_over_root_b = 0.0
            self.mean = 0.0
            self.variance = 0.0
            self.err_variance = np.inf
            self.err_mean = np.inf
            return
        
        # select only events in the min interval range
        mask = np.logical_and(
            self.range[0] <= invmass, invmass <= self.range[1])

        # evaluate statistical properties of events in min interval range
        self.mean = np.average(invmass[mask], weights=weights[mask])
        self.variance = np.average(
            np.power((invmass[mask] - self.mean), 2), weights=weights[mask])
        self.err_mean = np.sqrt(self.variance)/np.sum(mask)
        self.err_variance = self.get_err_variance(invmass[mask], weights[mask])
        signal_weights = weights[np.logical_and(is_signal,mask)]
        bkg_weights = weights[np.logical_and(np.logical_not(is_signal),mask)]
        if len(bkg_weights) > 0:
            self.s_over_root_b = np.sum(signal_weights)/np.sqrt(np.sum(bkg_weights))
        else:
            self.s_over_root_b = 0

    def get_err_variance(self, x, w):
        """ uncertainty on variance is (m4 - m2^2) / (4 n m2)"""
        if self.err_mean == np.inf:
            return np.inf
        if sum(w) == 0:
            return np.inf

        m4 = np.average(np.power(x - self.mean, 4), weights=w)
        return (m4 - self.variance**2)/(4 * len(x) * self.variance)

    def weighted_quantile(self, mass, weights, quantiles):
        """ calculates unbinned interval containing target percent of events """

        if len(mass) == 0: return []
        quantiles = [0.5-(quantiles/2), 0.5+(quantiles/2)]
        weighted_quantiles = np.divide(np.subtract(
            np.cumsum(weights), 0.5 * weights), np.sum(weights))
        return np.interp(quantiles, weighted_quantiles, mass)
