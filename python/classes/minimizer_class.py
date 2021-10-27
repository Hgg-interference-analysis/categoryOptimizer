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
        self.min_mva = min(data['diphoton_mva'].values)
        self.max_mva = max(data['diphoton_mva'].values)
        self.cats = []
        self.boundaries = []
        self.bounds = []
        self.minimum = np.inf
        self.optimal_boundaries = []

    def update_bounds(self):
        """ updates the rectangular bounds on the boundary locations """
        self.bounds = []
        for i in range(self.num_cats):
            self.bounds.append((self.min_mva + 0.05*float(i), self.max_mva - 0.05*float(self.num_cats-i)))

    def create_categories(self, use_bounds = False):
        """ creates categories for the minimizer """
        if len(self.cats) > 0:
            self.cats = []

        if not use_bounds:
            self.boundaries = []
            self.boundaries.append(self.min_mva)
            for i in range(self.num_cats-1):
                self.boundaries.append(rand.uniform(self.bounds[i+1][0], self.bounds[i+1][1]))

            self.boundaries.sort()
            for i in range(self.num_cats):
                lower, upper = self.boundaries[i], self.boundaries[i+1] if i != self.num_cats - 1 else 1.0
                invmass = self.data[self.data['diphoton_mva'].between(lower, upper)]
                self.cats.append(mva_category(np.array(invmass['CMS_hgg_mass'].values), np.array(invmass['weight'].values), lower, upper))

        else:
            for i in range(len(self.boundaries)):
                lower,upper = self.boundaries[i], self.boundaries[i+1] if i+1 < len(self.boundaries) else 1.0
                invmass = self.data[self.data['diphoton_mva'].between(lower, upper)]
                self.cats.append(mva_category(np.array(invmass['CMS_hgg_mass'].values), np.array(invmass['weight'].values), lower, upper))
        

    def target(self, boundaries):
        """ target function to minimize """
        self.boundaries = boundaries
        self.boundaries.sort()
        self.create_categories(use_bounds=True)
        total_resolution = 0
        for cat in self.cats:
            total_resolution += cat.resolution
        return np.sqrt(total_resolution)

    def optimize_boundaries(self):
        """ optimizes the location of the boundaries by minimizing the target function """
        eps = 0.001
        optimum = minimize(self.target,
                           self.boundaries,
                           method='L-BFGS-B',
                           bounds=self.bounds,
                           options={'eps':eps})
        return optimum.fun, optimum.x
        
    def run(self):
        """ runs the minimizer """
        for i in range(self.num_tests): 
            print(f'iter {i}')
            self.update_bounds()
            self.create_categories(use_bounds=False)
            print(self.bounds)
            val, optimum = self.optimize_boundaries()
            print(val, optimum)
            if val < self.minimum:
                self.optimal_boundaries = optimum.copy()
                self.minimum=val

