#!/usr/bin/env python3

"""
This is the control script for the diphoton mva boundary optimizer

Data should be taken from output from a dumper in flashgg.
"""

import argparse as ap
from numpy import log
import pandas as pd
import sys

import python.helpers.helper_optimize as helper_optimize
from python.classes.minimizer_class import minimizer
import python.plotters.plot as plot
import logging


def main():

    # option management
    parser = ap.ArgumentParser(description="Diphoton MVA Boundary Optimizer")

    parser.add_argument("-i", "--inputFile", required=True,
                        help='input config file containing dataset and tree whose resolution will be the target of the minimization')
    parser.add_argument("-o", "--output",
                        help='output tag used to generate output file names')
    parser.add_argument("-n", "--num-bounds", default=4, type=int, dest='n_bounds',
                        help='number of categories to optimize')
    parser.add_argument("--num-iter", default=10, type=int, dest='iters',
                        help='number of times to run the minimizer')
    parser.add_argument("--lumi-scale", default=[1.,1.,1.], type=float, dest='lumi', nargs="*", 
                        help="luminosity of this dataset, to be applied to signal")
    parser.add_argument("--lumi-scale-bkg", default=[1.,1.,1.], type=float, dest='lumi_bkg', nargs="*", 
                        help="luminosity of this dataset, to be applied to background")
    parser.add_argument("--bkg-scale", default=1., type=float, dest='bkg_scale',
                        help="scale for background file weights, use python/utils/derive_bkg_scale.py")
    parser.add_argument("--xcheck", default=False, action='store_true',
                        help='plots variables for each file')
    parser.add_argument("--plot", default=False, action='store_true',
                        help='makes the stack plots')
    parser.add_argument("-b", "--boundaries", type=float, default=[], nargs="*", 
                        help="Adds boundaries to the stack plot")
    parser.add_argument("--log", type=str, default=None, required=True,
                        help="name for log file")
    args = parser.parse_args()

    logging.basicConfig(filename=args.log, level='INFO')
    helper_optimize.print_setup(sys.argv, args.inputFile)

    df_sig, df_bkg = helper_optimize.extract_data(args)

    if args.xcheck or args.plot:
        return

    df_data = pd.concat([df_sig, df_bkg])
    # hand the data off the minimizer
    my_minimizer = minimizer(df_data, args.n_bounds, num_tests=args.iters, seed=args.boundaries)
    my_minimizer.run()

    opt_bounds = [round(x, 6) for x in my_minimizer.optimal_boundaries]
    opt_res = 100*round(my_minimizer.res[my_minimizer.minimum], 6)
    opt_res_err = round(100*my_minimizer.min_unc, 4)
    opt_sorb = my_minimizer.s_over_root_b
    log_string = f'the optimal boundaries are {opt_bounds}'
    logging.info(log_string)
    log_string = f'The associated resolution is {opt_res} +/- {opt_res_err}'
    logging.info(log_string)
    log_string = f'The s.o.r.b. is {round(opt_sorb,3)}'
    logging.info(log_string)


if __name__ == "__main__":
    main()
    