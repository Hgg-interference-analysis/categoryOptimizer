"""
Class for the minimizer 
"""
from matplotlib.pyplot import xcorr
import numpy as np
import random as rand
from scipy.optimize import minimize
import logging

from python.classes.category_class import mva_category


class minimizer:

    def __init__(self, data, num_cats, num_tests=10, seed=None) -> None:
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
        self.min_mva = min(self.mva)+0.025
        self.max_mva = max(self.mva)
        self.cats = []
        self.boundaries = []
        self.bounds = []
        self.minimum = np.inf
        self.s_over_root_b = 0
        self.optimal_boundaries = []
        self.res = {}
        self.res_uncs = {}
        self.signal_strengths = {}
        self.seed = seed
        self.first = True
        self.sort_data()

    def get_combined_resolution(self):
        # returns the combined relative resolution
        sigma = 0
        sigma_weight = 0
        mean = 0
        mean_weight = 0
        for cat in self.cats:
            sigma += cat.variance / cat.err_variance
            sigma_weight += 1 / cat.err_variance
            mean += cat.mean / cat.err_mean**2
            mean_weight += 1 / cat.err_mean**2
        if mean_weight == 0:
            return 999, 999
        w_mean = mean/mean_weight
        w_mean_unc = np.sqrt(1/mean_weight)
        w_sigma = np.sqrt(sigma/sigma_weight)
        w_sigma_unc = np.sqrt(1/sigma_weight)
        res = w_sigma/w_mean
        res_err = res * np.sqrt((w_mean_unc/w_mean)**2 
                                + (w_sigma_unc/w_sigma)**2)
        
        return res, res_err

    def get_sorb(self):
        # returns the quad added S over root B for each category
        return np.sum(cat.sig**2/cat.bkg for cat in self.cats)

    def sort_data(self):
        """ sorts data for faster quantile evaluation """
        sorter = np.argsort(self.mass)
        self.mass = self.mass[sorter]
        self.weights = self.weights[sorter]
        self.mva = self.mva[sorter]
        self.sig_bkg = self.sig_bkg[sorter]

    def update_bounds(self):
        """ updates the rectangular bounds on the boundary locations """
        self.bounds = []
        for i in range(self.num_cats):
            self.bounds.append((self.min_mva + 0.1*float(i),
                               self.max_mva - 1/float((3*i)**(3*i)+2)))

    def create_categories(self, use_bounds=False):
        """ creates categories for the minimizer """
        if len(self.cats) > 0:
            self.cats = []

        if not use_bounds:
            self.boundaries = []
            self.boundaries = [rand.uniform(
                self.bounds[i][0], self.bounds[i][1]) for i in range(self.num_cats)]
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
                self.cats.append(
                    mva_category(self.mass[mask], 
                                 self.weights[mask], 
                                 self.sig_bkg[mask], 
                                 )
                                 )

    def target(self, boundaries):
        """ target function to minimize """

        # don't allow boundaries that are too small
        diffs = [abs(boundaries[i]-boundaries[i+1])
                 for i in range(len(boundaries)-1)]
        if any(x <= 0.01 for x in diffs):
            self.res[999] = 999
            self.res_uncs[999] = 999
            self.signal_strengths[999] = 999
            return 999

        # assign and sort the boundaries, and create the categories
        self.boundaries = boundaries 
        self.boundaries.sort()
        self.create_categories(use_bounds=True)
            
        res, res_err = self.get_combined_resolution()
        sorb = self.get_sorb()
        ret = 999
        if sorb != 0:
            ret = 1000*res/sorb
        self.res[ret] = res
        self.res_uncs[ret] = res_err
        self.signal_strengths[ret] = sorb


        #print(f"{[round(x,4) for x in boundaries]} {round(ret,4)} {[round(100*np.sqrt(cat.variance)/cat.mean,4) for cat in self.cats]} {[round(cat.sig/np.sqrt(cat.bkg),2) for cat in self.cats]}")
        if not ret:
            self.max_mva = self.max_mva-0.001
            self.bounds[-1][1] = self.max_mva

        return ret

    def optimize_boundaries(self):
        """ optimizes the location of the boundaries by minimizing the target function """
        eps = 0.001
        optimum = minimize(self.target,
                           self.boundaries,
                           method='L-BFGS-B',
                           bounds=self.bounds,
                           options={'eps': eps}
                           )
        return optimum.fun, self.res[optimum.fun], self.res_uncs[optimum.fun], optimum.x, self.signal_strengths[optimum.fun]

    def run(self):
        """ runs the minimizer """
        for i in range(self.num_tests):
            self.update_bounds()
            self.res_uncs = {}
            if self.seed and i==0:
                self.boundaries = self.seed
                self.create_categories(use_bounds=True)
            else:
                self.create_categories(use_bounds=False)
            fun, val, val_unc, optimum, s_over_root_b = self.optimize_boundaries()
            if fun < self.minimum:
                optimum.sort()
                logging.info(f' iter{i}: {100*round(val,6)} +/- {100*round(val_unc,6)}, {round(s_over_root_b,3)}, {optimum}')
                #logging.info(f'{[cat.sig/np.sqrt(cat.bkg) for cat in self.cats]}')
                self.optimal_boundaries = optimum.copy()
                self.minimum = fun
                self.min_unc = val_unc
                self.s_over_root_b = s_over_root_b
