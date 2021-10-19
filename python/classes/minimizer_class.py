
import random as rand
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

    def create_categories(self):
        for i in range(num_cats):
            lower = 0
            upper = 0.1
            if len(self.cats) == 0:
                lower = rand.uniform(self.min_mva, self.max_mva - 0.5)
            else:
                lower = rand.uniform(self.cats[-1].upper, self.max_mva)
            
            upper = rand.uniform(lower, self.max_mva)

            self.cats.append(mva_category(lower,upper))
        
        self.cats[-1].set_upper_bound(1.0)

    def optimize_boundaries(self):
        
        
    def run(self):
        """ starts the minimizer """
        self.create_categories()
        self.optimize_boundaries()

