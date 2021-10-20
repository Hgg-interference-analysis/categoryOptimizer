#!/usr/bin/env python3

"""
This is the control script for the diphoton mva boundary optimizer

Data should be taken from output from a dumper in flashgg.

Usage/Order of Operations:

"""

import argparse as ap
import pandas as pd

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

    df_data = helper_optimize.extract_data(args.input_file)

    if args.xcheck:
        plot.plot(df_data)
        return

    #hand the data off the minimizer
    my_minimizer = minimizer(df_data, args.num_cats, num_test=args.num_tests, method=args.method)
    my_minimizer.run()

    print(f'the optimal boundaries are {my_minimizer.optimal_boundaries}')


if __name__ == "__main__":
    main()
