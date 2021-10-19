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
    keep_cols = ['DiphotonMVA','CMS_hgg_mass','leadEta','subleadEta','leadPt','subleadPt','lead_R9','sublead_R9','leadIDMVA','subleadIDMVA','decorrSigmaM','weight']
    keep_cols_bkg = ['DiphotonMVA','CMS_hgg_mass','leadEta','subleadEta','leadPt','subleadPt','lead_R9','sublead_R9','leadIDMVA','subleadIDMVA','decorrSigmaM','weight','gen_lead_pt','gen_sublead_pt']
    drop_cols = ['leadEta','subleadEta','leadPt','subleadPt','leadIDMVA','subleadIDMVA']
    

    #open and load the data
    config = open(args.inputFile,'r').readlines()
    config = [x.strip() for x in config]

    files = []

    for line in config:
        line_list = line.split('\t')
        print("[INFO] opening {} as dataframe".format(line_list[1]))
        df = pd.DataFrame()
        if line_list[0].find('Back') != -1 or line_list[0].find('back') != -1:
            df = up.open(line_list[1])[line_list[0]].pandas.df(keep_cols_bkg)
        else:
            df = up.open(line_list[1])[line_list[0]].pandas.df(keep_cols)
        #preselection
        files.append(df)

    if args.xcheck:
        plot.plot(files)
        return
        
    data = pd.concat(files)
    data.drop(drop_cols,axis=1,inplace=True)

    #hand the data off the minimizer
    bound_opt.minimize_target(data, args.num_bounds, args.num_iter,args.method)
    

main()
