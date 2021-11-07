import numpy as np
from numpy.core.shape_base import stack
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import iqr

def color_gradient(n):       
    """ returns a list of rgb values scaled from blue to green to red """    

    colors = []                                                                 
    for i in range(n):                                                          
        blue = (1 - 2*float(i)/float(n-1))                                      
        if blue < 0:                                                            
            blue = 0.                                                           
        red = float(i)/float(n-1) if blue == 0 else 0.                          
        green = 0.0                                                             
        if i == (n-1)/2:                                                        
            green = 1.                                                          
        elif i < (n-1)/2:                                                       
            green = 2 * float(i)/float(n-1)                                     
        else:                                                                   
            green = 2 - 2*float(i)/float(n-1)                                   
                                                                                
        colors.append((red,green,blue))                                         
                                                                                
    return colors

def get_bin_uncertainties(bins, values, weights):
    """ calculates the uncertainty of weighted bins """

    ret = []
    for i in range(len(bins)-1):
        val_mask = np.logical_and( bins[i] <= values, values < bins[i+1])
        ret.append(np.sqrt(np.sum(np.power(weights[val_mask], 2))))

    return np.array(ret)