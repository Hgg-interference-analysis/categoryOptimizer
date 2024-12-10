"""
Class for the minimizer 
"""
from matplotlib.pyplot import xcorr
import numpy as np
import random as rand
from scipy.optimize import minimize
import logging
import itertools

from python.classes.category_class_new import pt_category


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
        self.min_pt = 15.0
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


    def find_best_category(self):
        """ find the category with the smallest value of the loss function """

        # if there are already categories, remove them
        if len(self.cats) > 0:
            self.cats = []
   
        #  create boundaries with a specific step size:
        print("creating boundaries...")
        self.boundaries = []
        boundary = self.min_pt
        #while boundary < self.max_pt:
        while boundary < 130:  # very less events after ~120 GeV diphoton pT -- no point in creating categories
            self.boundaries.append(boundary)
            boundary += 5.0
        #print(f'all possible boundaries: {self.boundaries}')
        # create all possible combinations of boundaries
        partitions = itertools.combinations(self.boundaries, self.num_cats)
        final_loss = 9999.
        print("creating categories and calculating loss for each category")
        for partition in partitions:
            #print(f'partition: {partition}')
            invalid_flag = False
            self.cats = [] # reset categories
            for i in range(self.num_cats):
                lower = partition[i]
                upper = partition[i+1] if i+1 < self.num_cats else self.max_pt
                #print(f'lower boundary: {lower}')
                #print(f'upper boundary: {upper}')
                mask = np.logical_and(lower <= self.pt, self.pt < upper)
                self.cats.append(pt_category(lower, upper, self.mass[mask], 
                                self.weights[mask], 
                                self.sig_bkg[mask], 
                    )
                )
                if self.cats[-1].invalid:
                    invalid_flag = True
                    break
                    
            if invalid_flag:
                continue

            loss = self.calc_loss()
            #print(f'loss: {loss}')
            if loss < final_loss:
                final_loss = loss
                optimum_partition = partition
                final_cats = self.cats
                #print(f'optimum_partition: {optimum_partition}')
        
        return final_loss, optimum_partition, final_cats


    def calc_loss(self):
        """ target function to minimize (for diphoton pT boundary optimization)"""
        #self.boundaries = boundaries 
        #self.create_categories(use_bounds=True)
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
            #print(f'res: {res}')
            res_err = res * np.sqrt((mean_err/mean)**2 
                                + (sigma_err/sigma)**2)
            
            # get significance
            s_o_r_b = cat.sorb
            #print(f'sorb: {s_o_r_b}')
            lfunc += (res / s_o_r_b)
  
        #print(f'lfunc: {lfunc}')
        ret = 999
        if lfunc != 0 and not np.isnan(lfunc):
            ret = 1000*lfunc   
        #print(f'ret: {ret}')
        return ret


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
        print("running the minimizer")
        optimum_loss, optimum_bounds, optimum_cats = self.find_best_category()

        #if (fun < self.minimum) and fun != 999:
            #optimum.sort()
            #logging.info(f' iter{i}: {100*round(val,6)} +/- {100*round(val_unc,6)}, {round(s_over_root_b,3)}, {optimum}, {fun}')
        self.optimal_boundaries = optimum_bounds
        self.minimum = optimum_loss
        #logging.info(f'{optimum,4}, {fun}')
        print ("minimum = ", self.minimum)
        for j,cat in enumerate(optimum_cats):
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
                    

