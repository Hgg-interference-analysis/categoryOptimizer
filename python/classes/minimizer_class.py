"""
Class for the minimizer 
"""
import random as rand
import numpy as np

from scipy.optimize import minimize

from python.classes.category_class import mva_category

class minimizer:

    def __init__(self, data, num_cats, num_tests=10, method='series') -> None:
        self.data = data 
        self.num_cats = num_cats
        self.num_tests = num_tests
        self.method = method
        self.min_mva = min(data['DiphotonMVA'].values)
        self.max_mva = 1.0
        self.cats = []
        self.boundaries = []
        self.bounds = []
        self.minimum = np.inf
        self.optimal_boundaries = []

    def update_bounds(self):
        """ updates the rectangular bounds on the boundary locations """
        for i in range(len(self.cats)):
            if i == 0:
                self.bounds.append((self.min_mva, self.cats[i+1].lower))
            elif 0 < i < len(self.cats)-1:
                self.bounds.append((self.cats[i-1].upper, self.cats[i+1].lower))
            else:
                self.bounds.append((self.cats[i-1].upper, 1.0))

    def create_categories(self, use_bounds = False):
        """ creates categories for the minimizer """
        if len(self.cats) > 0:
            self.cats = []

        if not use_bounds:
            for i in range(num_cats):
                lower = 0.
                upper = 0.1
                if len(self.cats) == 0:
                    lower = rand.uniform(self.min_mva, self.max_mva - 0.5)
                    self.boundaries.append(lower)
                else:
                    lower = rand.uniform(self.cats[-1].upper, self.max_mva)
            
                upper = rand.uniform(lower, self.max_mva) if i != num_cats-1 else 1.0
                invmass = self.data[self.data['DiphotonMVA'].between(lower, upper)]
                self.cats.append(mva_category(np.array(invmass['CMS_hgg_mass'].values), np.array(invmass['weight'].values), lower, upper))
                if i != num_cats-1:
                    self.boundaries.append(upper)

        else:
            for i in range(len(self.boundaries)):
                lower = self.boundaries[i]
                upper = self.boundaries[i+1] if i+1 < len(self.boundaries) else 1.0
                invmass = self.data[self.data['DiphotonMVA'].between(lower, upper)]
                self.cats.append(mva_category(np.array(invmass['CMS_hgg_mass'].values), np.array(invmass['weight'].values), lower, upper))
        

    def target(self):
        """ target function to minimize """
        self.create_categories(use_bounds=True)
        self.update_bounds()
        total_resolution = 0
        for cat in self.cats:
            total_resolution += cat.resolution
        
        return np.sqrt(total_resolution)

    def optimize_boundaries(self):
        """ optimizes the location of the boundaries by minimizing the target function """
        eps = 0.0001
        optimum = minimize(self.target,
                           guess=self.boundaries
                           method='L-BFGS-B',
                           bounds=self.bounds,
                           options={'eps':eps})
        return optimum.fun, optimum.x
        
    def run(self):
        """ runs the minimizer """
        while self.num_test > 0:
            self.num_tests -= 1
            self.create_categories(use_bounds=False)
            self.update_bounds()
            val, optimum = self.optimize_boundaries()
            if val < self.minimum:
                self.optimal_boundaries = optimum.copy()

