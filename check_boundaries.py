#!/usr/bin/env python3

import argparse as ap
import pandas as pd
import numpy as np
import math

import python.helpers.helper_optimize as helper_optimize
from python.classes.minimizer_class import minimizer
import python.plotters.plot as plot


def transform(b):
    """ transforms the boundaries into the new phase space """

    for i,val in enumerate(b):
        b[i] = 1 / (1 + math.exp(0.5 * math.log(2./ (val + 1) - 1)))

    return b

    

def main():

    # option management
    parser = ap.ArgumentParser(description="Diphoton MVA Boundary Checker")

    parser.add_argument("-i", "--inputFile", required=True,
                        help='input config file containing dataset and tree whose resolution will be the target of the minimization')
    parser.add_argument("-b", "--boundaries", type=float, default=[], nargs="*", 
                        help='boundaries to check')
    parser.add_argument("--transform", default=False, action="store_true",
                        help="will transform the values into the new mva space")
    parser.add_argument("--lumi-scale", default=1., type=float, dest='lumi',    
                        help="luminosity of this dataset")
    parser.add_argument("--bkg-scale", default=1., type=float, dest='bkg_scale',
                        help="scale for background file weights, use python/utils/derive_bkg_scale.py")

    args = parser.parse_args()

    old_bounds = args.boundaries.copy()
    old_bounds.sort()
    if args.transform:
        args.boundaries = transform(args.boundaries)

    args.boundaries.sort()
    args.xcheck = False
    args.plot = False
    args.output = '' 

    df_sig, df_bkg = helper_optimize.extract_data(args)

    df_data = pd.concat([df_sig, df_bkg])
    # hand the data off the minimizer
    my_minimizer = minimizer(df_data, len(args.boundaries))

    my_minimizer.boundaries = args.boundaries
    my_minimizer.update_bounds()
    my_minimizer.create_categories(use_bounds=True)

    res = my_minimizer.target(my_minimizer.boundaries)
    res_err = my_minimizer.res_uncs[res]
    sorb = my_minimizer.signal_strengths[res]

    print()
    print(f'The boundaries we are testing are {[round(x,4) for x in old_bounds]}')
    if args.transform:
        print(f'The transformed values are {[round(x,4) for x in args.boundaries]}')
    print(f'The associated resolution is {round(100*res,4)} +/- {round(100*res_err,4)} %')
    print(f'The s.o.r.b. is {round(sorb,2)}')


if __name__ == "__main__":
    main()