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
        #self.mva = np.array(data['diphoton_transformedMva'].values)  
        self.mva = np.concatenate([data['diphoton_mva'].dropna().values, data['diphoton_transformedMva'].dropna().values])
        self.pt = np.array(data['diphoton_pt'].values)                    
        self.mass = np.array(data['CMS_hgg_mass'].values)
        self.weights = np.array(data['weight'].values)
        self.sig_bkg = np.array(data['is_signal'].values)
        mass_mask = np.logical_and(115 < self.mass, self.mass < 135)  
        mva_mask =  self.mva > 0.6 
        self.mass = self.mass[mass_mask & mva_mask]
        self.mva = self.mva[mass_mask & mva_mask]
        self.pt = self.pt[mass_mask & mva_mask]
        self.weights = self.weights[mass_mask & mva_mask]
        self.sig_bkg = self.sig_bkg[mass_mask & mva_mask]
        self.num_cats = num_cats
        self.num_tests = num_tests
        self.min_mva = min(self.mva)+0.025    
        self.max_mva = max(self.mva)-0.025  
        self.min_pt = 15
        self.max_pt = max(self.pt)  
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
        self.count = 0
    
    def get_combined_resolution(self):
        """ returns the combined resolution """
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
        """ returns the quad-added sorb for each category"""
    
        return np.sum(cat.sig**2/cat.bkg for cat in self.cats)

    def sort_data(self):
        """ sorts data for faster quantile evaluation """

        sorter = np.argsort(self.mass)
        self.mass = self.mass[sorter]
        self.weights = self.weights[sorter]
        self.mva = self.mva[sorter]
        self.sig_bkg = self.sig_bkg[sorter]
        self.pt = self.pt[sorter]

    def update_bounds(self):
        """ updates the rectangular bounds on the boundary locations """

        self.bounds = []
        for i in range(len(self.boundaries)):
            b = (self.min_mva if i == 0 else self.boundaries[i-1]+0.001, self.max_mva if i == len(self.boundaries)-1 else self.boundaries[i+1]-0.001)
            self.bounds.append(b)
        print("Updated bounds :", self.bounds)


    def create_categories(self, use_bounds=False):
        """ creates categories for the minimizer """

        # if there are already categories, remove them
        if len(self.cats) > 0:
            self.cats = []
   
        # create the categories
        for i in range(len(self.boundaries)):
            if len(self.boundaries) != self.num_cats:
                print("Number of categories do not match the number of boundaries given. Creating only {} categories".format(len(self.boundaries)))
            lower = self.boundaries[i]
            upper = self.boundaries[i+1] if i+1 < len(self.boundaries) else self.max_pt
            mask = np.logical_and(lower <= self.pt, self.pt < upper)
            self.cats.append(mva_category(lower, upper, self.mass[mask], 
                                self.weights[mask], 
                                self.sig_bkg[mask], 
                    )
                )


    def target_for_pt_opt(self, boundaries):
        """ target function to minimize (for diphoton pT boundary optimization)"""
        self.boundaries = boundaries 
        self.create_categories(use_bounds=True)
        lfunc = 0

        """ relative mass resolution """
        for cat in self.cats:
            sigma = np.sqrt(cat.variance)
            sigma_err = np.sqrt(cat.err_variance)    
            mean = cat.mean
            mean_err = cat.err_mean
            mean_weight = 1 / cat.err_mean**2
            if mean_weight == 0:
                return 999
            
            res = sigma/mean
            res_err = res * np.sqrt((mean_err/mean)**2 
                                + (sigma_err/sigma)**2)
            
            # get significance
            s_o_r_b = cat.sorb
            lfunc += (res / s_o_r_b)
  
        #print(lfunc)
        ret = 999
        if lfunc != 0 and not np.isnan(lfunc):
            ret = 100*lfunc   
        
        #print(ret)
        return ret

    #def callback_func(x):
    #    if self.count_nan >= 10:
    #        return True

    
    def target(self, boundaries):
        """ target function to minimize (for diphoton MVA score boundary optimization)"""
        # assign and sort the boundaries, and create the categories
        self.boundaries = boundaries 
        self.create_categories(use_bounds=True)
        
            
        # evaluate properties of the current boundaries
        res, res_err = self.get_combined_resolution()
        sorb = self.get_sorb()
        ret = 999
        if sorb != 0 and not np.isnan(sorb):
            ret = 1000*res/sorb    
        self.res[ret] = res
        self.res_uncs[ret] = res_err
        self.signal_strengths[ret] = sorb

        #print(f"{[round(x,6) for x in boundaries]} {round(ret,4)} {[round(100*np.sqrt(cat.variance)/cat.mean,4) for cat in self.cats]} {[round(cat.sig/np.sqrt(cat.bkg),2) for cat in self.cats]}")

        return ret


    def optimize_boundaries(self):
        """ optimizes the location of the boundaries by minimizing the target function """
        optimum = minimize(self.target_for_pt_opt,
                           self.boundaries,
                           method='Nelder-Mead', # uses a greedy algorithm, it's good to run this many times to find the minimum
                           bounds=self.bounds, 
                           tol=1e-4, 
                           options={'fatol': 1e-4} 
                           )
        
        if np.isnan(optimum.fun):
            #return np.inf, 999, 999, optimum.x, [-1 for x in range(4)]
            return 999, optimum.x
        #return optimum.fun, self.res[optimum.fun], self.res_uncs[optimum.fun], optimum.x, self.signal_strengths[optimum.fun]
        return optimum.fun, optimum.x


    def run(self):
        """ runs the minimizer """

        # loop over number of iterations
        for i in range(self.num_tests):

            self.res_uncs = {}

            # if this is the first run, and a seed is provided, use the seed boundaries
            if self.seed and i==0:
                self.boundaries = self.seed
            else:
                self.boundaries = []
                self.boundaries = [round(rand.uniform(self.min_pt, 110),0) for i in range(self.num_cats)]   #restricted the initial guess to not go over 110 \
                                                                                                   #(very less events and no intereference effects after that).
                self.boundaries.sort()
            invalid_ = False
            #print("Boundaries :", self.boundaries)
            for k in range(len(self.boundaries) - 1):
                if self.boundaries[k+1] - self.boundaries[k] < 10:
                    invalid_ = True
            if invalid_:
                self.count+=1
                continue 
            # update the bounds based on the current boundaries
            #self.update_bounds()
            self.bounds = []
            for n in range(self.num_cats):
                self.bounds.append((15,115))      ## fixed bounds for each boundary
                
            # optimize the current boundaries
            #fun, val, val_unc, optimum, s_over_root_b = self.optimize_boundaries()
            fun, optimum = self.optimize_boundaries()

            # if the current boundaries outperformed the last boundaries, make a note of it
            if (fun < self.minimum) and fun != 999:
                #optimum.sort()
                #logging.info(f' iter{i}: {100*round(val,6)} +/- {100*round(val_unc,6)}, {round(s_over_root_b,3)}, {optimum}, {fun}')
                self.optimal_boundaries = optimum.copy()
                self.minimum = fun
                #logging.info(f'{optimum,4}, {fun}')
                print ("minimum = ", self.minimum)
                logging.info(' iter{}:'.format(i))
                for j,cat in enumerate(self.cats):
                    print("cat {}".format(j))
                    print("s.o.r.b = ", round(cat.sorb, 4) )
                    print("lower boundary = ", round(cat.lower_boundary, 4))
                    print("upper boundary = ", round(cat.upper_boundary, 4))
                    print("total no. of events = ", round(cat.nEvents, 1))
                    print("no. of signal events = ", round(cat.nSig, 1))
                    print("no. of bkg events = ", round(cat.nBkg, 1))
                    print("sig statistics = ", round(cat.nSigEvnts, 1))
                    print("bkg statistics = ", round(cat.nBkgEvnts, 1))
                    resolution = np.sqrt(cat.variance) / cat.mean
                    res_unc = resolution * np.sqrt((cat.err_mean/cat.mean)**2 
                                + (np.sqrt(cat.err_variance)/np.sqrt(cat.variance))**2)
                    print("resolution = ", round(resolution, 6))
                    print("error in resolution = ", round(res_unc, 6))
                
                    logging.info('cat{}: {} +/- {}, {}'.format(j, round(resolution, 4), round(res_unc, 6), round(cat.sorb, 3)))
                
                    

