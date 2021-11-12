
from matplotlib.pyplot import xcorr
import numpy as np
from scipy.stats import moment


class mva_category:
    """ class to describe diphoton mva categories """

    def __init__(self, invmass, weights, is_signal) -> None:
        self.range = self.weighted_quantile(invmass, weights, 0.683)
        mask = np.logical_and(
            self.range[0] <= invmass, invmass <= self.range[1])
        self.mean = np.average(invmass[mask], weights=weights[mask])
        self.variance = np.average(
            np.power((invmass[mask] - self.mean), 2), weights=weights[mask])
        # print(self.range, np.sqrt(self.variance), self.mean, np.sqrt(self.variance)/self.mean)
        self.err_variance = self.get_err_variance(invmass[mask], weights[mask])
        self.err_mean = np.sqrt(self.variance)/np.sum(mask)
        signal_weights = weights[np.logical_and(is_signal,mask)]
        bkg_weights = weights[np.logical_and(np.logical_not(is_signal),mask)]
        self.s_over_root_b = np.sum(signal_weights)/np.sqrt(np.sum(bkg_weights))

    def get_err_variance(self, x, w):
        """ uncertainty on variance is (m4 - m2^2) / (4 n m2)"""

        m4 = np.average(np.power(x - self.mean, 4), weights=w)
        return (m4 - self.variance**2)/(4 * len(x) * self.variance)

    def weighted_quantile(self, mass, weights, quantiles):
        """ calculates unbinned interval containing target percent of events """

        quantiles = [0.5-(quantiles/2), 0.5+(quantiles/2)]
        weighted_quantiles = np.divide(np.subtract(
            np.cumsum(weights), 0.5 * weights), np.sum(weights))
        return np.interp(quantiles, weighted_quantiles, mass)
