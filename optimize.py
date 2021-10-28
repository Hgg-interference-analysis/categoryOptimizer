#!/usr/bin/env python3

"""
This is the control script for the diphoton mva boundary optimizer

Data should be taken from output from a dumper in flashgg.

Usage/Order of Operations:

"""

import argparse as ap
import pandas as pd
import sys

import python.helpers.helper_optimize as helper_optimize
from python.classes.minimizer_class import minimizer

###############################################################################
def main():

    #option management
    parser = ap.ArgumentParser(description="Diphoton MVA Boundary Optimizer")

    parser.add_argument("-i","--inputFile", required=True,
                        help='input config file containing dataset and tree whose resolution will be the target of the minimization')
    parser.add_argument("-o","--output",
                        help='output tag used to generate output file names')
    parser.add_argument("-n","--num-bounds", default=4, type=int, dest='n_bounds',
                        help='number of categories to optimize')
    parser.add_argument("--num-iter", default=10, type=int, dest='iters',
                        help='number of times to run the minimizer')
    parser.add_argument("--method", default='series', type=str,
                        help='method to use in minimization, options are "series" (default) or "parallel"')
    parser.add_argument("--xcheck", default=False, action='store_true',
                        help='plot a bunch of stuff')

    args = parser.parse_args()

    helper_optimize.print_setup(sys.argv, args.inputFile)

    df_sig, df_bkg = helper_optimize.extract_data(args.inputFile, args.xcheck)
    if args.xcheck:
        return

    df_data = pd.concat([df_sig, df_bkg])
    #hand the data off the minimizer
    my_minimizer = minimizer(df_data, args.n_bounds, num_tests=args.iters)
    my_minimizer.run()

    print(f'the optimal boundaries are {my_minimizer.optimal_boundaries}\nThe associated resolution is {my_minimizer.minimum}')


if __name__ == "__main__":
    main()
