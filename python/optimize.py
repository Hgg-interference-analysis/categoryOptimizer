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
                        

    args = parser.parse_args()

    cmd = ''
    for arg in sys.argv:
        if ' ' in arg:
            cmd += '"{}" '.format(arg)
        else:
            cmd +=  "{} ".format(arg)

    print("#"*40)
    print("[INFO] Welcome to the diphoton mva boundary optimizer")
    print("[INFO] the command you ran was: {}".format(cmd))
    print("[INFO] the specified config file is: {}".format(args.inputFile))

    #columns needed for minimzation are hardcoded
    keep_cols = ['DiphotonMVA','CMS_hgg_mass']
    

    #open and load the data
    config = open(args.inputFile,'r').readlines()
    config = [x.strip() for x in config]

    files = []

    for line in config:
        line_list = line.split('\t')
        print("[INFO] opening {} as dataframe".format(line_list[1]))
        df = up.open(line_list[1])[line_list[0]].pandas.df(keep_cols)
        mva_mask = df['DiphotonMVA'].between(-0.5,999)
        invmass_mask = df['CMS_hgg_mass'].between(105,140)
        df = df.loc[invmass_mask&mva_mask]
        files.append(df)

    data = pd.concat(files)

    #hand the data off the minimizer
    bound_opt.minimize_target(data, args.num_bounds, args.num_iter)
    

main()
