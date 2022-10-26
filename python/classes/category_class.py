
import matplotlib.pyplot as plt
import numpy as np

MIN_EVENTS = 100

class mva_category:
    """ class to describe diphoton mva categories """

    def __init__(self, invmass, weights, is_signal, do_sm=True) -> None:

        # evaluate range of 68.3% min interval
        
        if sum(weights) < MIN_EVENTS:
            self.set_invalid()
            return

        # select only events in the min interval range
        self.range = self.weighted_quantile(invmass[is_signal], weights[is_signal], 0.683)
        #print(self.range)
        if len(self.range) == 0 or sum(weights) == 0:
            self.set_invalid()
            return

        mask = np.logical_and(
            self.range[0] <= invmass, invmass <= self.range[1])

        #hist, bins = np.histogram(invmass, bins=20, range = self.range, normed=True, weights=weights)
        #mids = [(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]

        #plt.scatter(mids,hist,marker='.', color='b')
        #plt.scatter(mids, hist1, marker='o', color='r')
        #plt.show()

        if sum(mask) == 0 or sum(weights[mask]) == 0:
            self.set_invalid()
            return


        # sorb stuff
        signal_weights = weights[np.logical_and(is_signal, mask)]
        bkg_weights = weights[np.logical_and(np.logical_not(is_signal), mask)]
        if len(bkg_weights) > 0:
            self.sig = np.sum(signal_weights)
            self.bkg = np.sum(bkg_weights)
        else:
            self.set_invalid()
            return
            self.bkg = np.inf
            self.sig = 0

        if do_sm:

            # evaluate statistical properties of events in min interval range
            self.mean = np.average(invmass[mask], weights=weights[mask])
            self.variance = np.average(
                np.power((invmass[mask] - self.mean), 2), weights=weights[mask])
            self.err_mean = np.sqrt(self.variance)/np.sum(mask)
            self.err_variance = self.get_err_variance(
                invmass[mask], weights[mask])

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

        if len(mass) == 0:
            return []
        quantiles = [0.5-(quantiles/2), 0.5+(quantiles/2)]
        weighted_quantiles = np.divide(np.subtract(
            np.cumsum(weights), 0.5 * weights), np.sum(weights))
        return np.interp(quantiles, weighted_quantiles, mass)


    def set_invalid(self):
        self.mean = 0.0
        self.variance = np.inf
        self.err_variance = np.inf
        self.err_mean = np.inf
        self.sig = 0
        self.bkg = np.inf
