"""
Boundary Optimizer
    Optimizes the boundaries with the number of specified categories
"""

import numpy as np
import pandas as pd
import statistics as stats
from scipy.optimize import minimize as minz
from scipy import stats as sci_stat

def min_interval(x, interval):
    if len(x['CMS_hgg_mass']) < 1:
        print("empty")
        return x
    min_val = np.percentile(x['CMS_hgg_mass'].values,50-(interval/2))
    max_val = np.percentile(x['CMS_hgg_mass'].values,50+(interval/2))
    mask = x['CMS_hgg_mass'].between(min_val,max_val)
    return x[mask]

def target(x):
    """ 
    Target Function returning the combined uncertainty on the invariant mass
    """
    global _num_cats
    global _data_
    global _min_
    global _max_
    global _min_bound_
    global _max_bound_

    x.sort() #boundaries must be in ascending order

    ret = 0
    np.append(x,_max_)

    for i in range(_num_cats_-1):
        mask = _data_['DiphotonMVA'].between(x[i],x[i+1])
        inv_mass = _data_[mask]
        reduced_inv_mass = np.array(min_interval(inv_mass,_min_interval_)['CMS_hgg_mass'].values)
        if len(reduced_inv_mass) > 1:
            ret += pow(np.std(reduced_inv_mass)/len(reduced_inv_mass),2)
        else:
            print('[ERROR] no events in bin ({},{})'.format(x[i],x[i+1]))
            ret += 1

    return np.sqrt(ret)

def minimize_target(data, num_cats, iterations):
    """Handles the minimization"""
    print("#"*40)
    print("[INFO][boundary_optimizer] Beginning minimization")
    print("[INFO][boundary_optimizer] There are {} categories".format(num_cats))
    print("[INFO][boundary_optimizer] {} points will be tested".format(iterations))

    global _num_cats_
    global _data_
    global _min_
    global _max_
    global _min_bound_
    global _max_bound_
    global _min_interval_

    _data_ = data
    _num_cats_ = num_cats
    _min_ = min(data['DiphotonMVA'].values)
    _max_ = 1.0
    _min_bound_ = min(data['DiphotonMVA'].values)
    _max_bound_ = max(data['DiphotonMVA'].values-.1)
    _min_interval_ = 68.3

    bounds = [(_min_bound_+i*(0.1),_max_bound_-(num_cats-i)*0.1) for i in range(num_cats)]

    optima = []
    for i in range(iterations):
        guess = np.random.uniform(low=_min_bound_*np.ones(num_cats),
                            high=_max_bound_*np.ones(num_cats)).ravel().tolist()
        guess.sort()
        eps = 0.0001
        optima.append(minz(target, np.array(guess), method='L-BFGS-B', bounds=bounds, options={'eps':eps}))

    best = 999
    ret = optima[0]
    for opt in optima:
        if opt.fun < best:
            best = opt.fun
            ret = opt

    print(ret)
