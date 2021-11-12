"""
Class for the minimizer 
"""
from matplotlib.pyplot import xcorr
import numpy as np
import random as rand
from scipy.optimize import minimize

from python.classes.category_class import mva_category


class minimizer:

    def __init__(self, data, num_cats, num_tests=10) -> None:
        self.mva = np.array(data['diphoton_mva'].values)
        self.mass = np.array(data['CMS_hgg_mass'].values)
        self.weights = np.array(data['weight'].values)
        self.sig_bkg = np.array(data['is_signal'].values)
        mass_mask = np.logical_and(115 < self.mass, self.mass < 135)
        self.mass = self.mass[mass_mask]
        self.mva = self.mva[mass_mask]
        self.weights = self.weights[mass_mask]
        self.sig_bkg = self.sig_bkg[mass_mask]
        self.num_cats = num_cats
        self.num_tests = num_tests
        self.min_mva = min(self.mva)
        self.max_mva = max(self.mva)
        self.cats = []
        self.boundaries = []
        self.bounds = []
        self.minimum = np.inf
        self.optimal_boundaries = []
        self.res_uncs = {}
        self.signal_strengths = {}
        self.sort_data()

    def sort_data(self):
        sorter = np.argsort(self.mass)
        self.mass = self.mass[sorter]
        self.weights = self.weights[sorter]
        self.mva = self.mva[sorter]
        self.sig_bkg = self.sig_bkg[sorter]

    def update_bounds(self):
        """ updates the rectangular bounds on the boundary locations """
        self.bounds = []
        for i in range(self.num_cats):
            self.bounds.append((self.min_mva + 0.05*float(i),
                               self.max_mva - 0.01*float(self.num_cats-i)))

    def create_categories(self, use_bounds=False):
        """ creates categories for the minimizer """
        if len(self.cats) > 0:
            self.cats = []

        if not use_bounds:
            self.boundaries = []
            self.boundaries = [self.min_mva if i == 0 else rand.uniform(
                self.bounds[i][0]+0.05, self.bounds[i][1]) for i in range(self.num_cats)]
            self.boundaries.sort()

            for i in range(self.num_cats):
                lower = self.boundaries[i],
                upper = self.boundaries[i + 1] if i != self.num_cats - 1 else 1.0
                mask = np.logical_and(lower <= self.mva, self.mva <= upper)
                self.cats.append(mva_category(
                    self.mass[mask], self.weights[mask], self.sig_bkg[mask]))

        else:
            self.boundaries.sort()
            for i in range(len(self.boundaries)):
                lower = self.boundaries[i]
                upper = self.boundaries[i+1] if i+1 < len(self.boundaries) else 1.0
                mask = np.logical_and(lower <= self.mva, self.mva <= upper)
                self.cats.append(mva_category(
                    self.mass[mask], self.weights[mask], self.sig_bkg[mask]))

    def target(self, boundaries):
        """ target function to minimize """
        diffs = [abs(boundaries[i]-boundaries[i+1])
                 for i in range(len(boundaries)-1)]
        if any(x <= 0.01 for x in diffs):
            self.res_uncs[999] = 999
            self.signal_strengths[999] = 999
            return 999
        self.boundaries = boundaries
        self.boundaries.sort()
        self.create_categories(use_bounds=True)
        sigma = 0
        sigma_weight = 0
        mean = 0
        mean_weight = 0
        sorb = 0
        for cat in self.cats:
            sigma += cat.variance / cat.err_variance
            sigma_weight += 1 / cat.err_variance
            mean += cat.mean / cat.err_mean**2
            mean_weight += 1 / cat.err_mean**2
            sorb += cat.s_over_root_b**2
        w_mean = mean/mean_weight
        w_mean_unc = np.sqrt(1/mean_weight)
        w_sigma = np.sqrt(sigma/sigma_weight)
        w_sigma_unc = np.sqrt(1/sigma_weight)
        res = w_sigma/w_mean
        res_err = res * np.sqrt((w_mean_unc/w_mean) **
                                2 + (w_sigma_unc/w_sigma)**2)
        self.res_uncs[res] = res_err
        self.signal_strengths[res] = np.sqrt(sorb)
        return res

    def optimize_boundaries(self):
        """ optimizes the location of the boundaries by minimizing the target function """
        eps = 0.001
        optimum = minimize(self.target,
                           self.boundaries,
                           method='L-BFGS-B',
                           bounds=self.bounds,
                           options={'eps': eps})
        return optimum.fun, self.res_uncs[optimum.fun], optimum.x, self.signal_strengths[optimum.fun]

    def run(self):
        """ runs the minimizer """
        for i in range(self.num_tests):
            print(f'iter {i}')
            self.update_bounds()
            self.res_uncs = {}
            self.create_categories(use_bounds=False)
            val, val_unc, optimum, s_over_root_b = self.optimize_boundaries()
            print(f'{100*round(val,6)} +/- {100*round(val_unc,6)}', s_over_root_b, optimum)
            if val < self.minimum:
                self.optimal_boundaries = optimum.copy()
                self.minimum = val
                self.min_unc = val_unc
                self.s_over_root_b = s_over_root_b
