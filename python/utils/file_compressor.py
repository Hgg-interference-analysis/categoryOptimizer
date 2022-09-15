#!/usr/bin/env python3

from tkinter import W
import numpy as np
import pandas as pd
import argparse as ap
import uproot3 as up
import matplotlib.pyplot as plt



def main():

    parser = ap.ArgumentParser(description="options for bkg scale derivation")
    parser.add_argument("-f", "--file", default=None, type=str,
                        help="root file to compress")
    parser.add_argument("-t","--tree", default=None, type=str,
                        help="name of tree to compress")
    parser.add_argument("-o","--output", default=None, type=str,
                        help="name of file to write")
    parser.add_argument("--keep", default=None, type=str,
                        help="comma separated list of variable to keep")
    parser.add_argument("--signal", default=False, action="store_true",
                        help="indicates whether you are processing signal")

    args = parser.parse_args()

    #process keep cols
    keep_cols = args.keep.split(",")

    #open the file using uproot
    df = up.open(args.file)[args.tree].pandas.df(keep_cols)

    #create new file
    out = up.create(args.output)

    if not args.signal:
        df['weight'] = np.multiply(df['weight'].values, df['Norm_SFs'].values)
        df.drop("Norm_SFs", axis=1, inplace=True)

    df['diphoton_mva'] = df["diphoMVANew"].values
    df.drop("diphoMVANew", axis=1, inplace=True)

    df['diphoton_mva'] = 1. / ( 1. + np.exp( 0.5*np.log( (2./(df['diphoton_mva'].values+1.)) - 1 ) ) )

    my_dict = {}
    for col in df.columns:
        my_dict[col] = "float64"
    out[args.tree] = up.newtree(my_dict)

    my_dict = {}
    for col in df.columns:
        my_dict[col] = np.array(df[col].values)

    print(my_dict)
    out[args.tree].extend(my_dict)
    out[args.tree].show()

    out.close()


if __name__ == "__main__":
    main()
