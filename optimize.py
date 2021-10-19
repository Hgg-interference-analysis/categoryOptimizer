#!/usr/bin/env python3

"""
This is the control script for the diphoton mva boundary optimizer

Data should be taken from output from flashgg.

Usage/Order of Operations:

"""

import argparse as ap
import numpy as np
import os 
import pandas as pd
import sys
import uproot as up

import boundary_optimizer as bound_opt
import plot

import python.helpers.helper_optimize as helper_optimize

###############################################################################
def main():

    #option management
    parser = ap.ArgumentParser(description="Diphoton MVA Boundary Optimizer")

    parser.add_argument("-i","--inputFile", required=True,
                        help='input config file containing dataset and tree whose resolution will be the target of the minimization')
    parser.add_argument("-o","--output",
                        help='output tag used to generate output file names')
    parser.add_argument("-n","--num_bounds", default=4, type=int,
                        help='number of categories to optimize')
    parser.add_argument("--num_iter", default=10, type=int,
                        help='number of times to run the minimizer')
    parser.add_argument("--method", default='series', type=str,
                        help='method to use in minimization, options are "series" (default) or "parallel"')
    parser.add_argument("--xcheck", default=False, action='store_true',
                        help='plot a bunch of stuff')

    args = parser.parse_args()

    helper_optimize.print_setup(sys.argv, args.inputFile)

    df_collection_data = helper_optimize.extract_data(args.input_file)

    if args.xcheck:
        plot.plot(df_collection_data)
        return
        
    data = pd.concat(df_collection_data)
    data.drop(drop_cols,axis=1,inplace=True)

    #hand the data off the minimizer
    bound_opt.minimize_target(data, args.num_bounds, args.num_iter,args.method)
    

if __name__ == "__main__":
    main()
