#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse as ap
import uproot3 as up
import matplotlib.pyplot as plt


def load_cfg_to_df(cfg_file):
    """ loads the cfg file into dataframes """

    config = open(cfg_file, 'r').readlines()
    config = [x.strip() for x in config]

    keep_cols = ['CMS_hgg_mass', 'diphoton_mva', 'weight']

    bkg_files = []
    bkg_titles = []
    data_files = []
    data_titles = []

    for line in config:
        line_list = line.split('\t')
        print("[INFO] opening {} as dataframe".format(line_list[1]))
        df = pd.DataFrame()
        if line_list[0].find('bkg') != -1:
            df = up.open(line_list[2])[line_list[1]].pandas.df(keep_cols)
            bkg_files.append(df)
            bkg_titles.append(line_list[3])
        else:
            df = up.open(line_list[2])[line_list[1]].pandas.df(keep_cols)
            data_files.append(df)
            data_titles.append(line_list[3])

    df_bkg = pd.concat(bkg_files)
    df_data = pd.concat(data_files)

    return df_data, df_bkg


def get_side_bands_hist(df, var, hrange=[100, 180], num_bins=160, _kIsData=False, lumi_scale=41.5):
    """ returns a numpy array of the sideband regions """

    mass = np.array(df[var].values)
    weights = lumi_scale * \
        np.array(df['weight'].values) if not _kIsData else np.ones(len(mass))

    side_band_low = 115
    side_band_high = 135
    mask = np.logical_or(mass < side_band_low, side_band_high < mass)
    mass = mass[mask]
    weights = weights[mask]

    this_var = (df[mask])[var]
    return np.histogram(this_var, bins=num_bins, range=hrange, weights=weights)


def plot(bins, data, bkg, xlabel, ylabel, title):
    """ plots data and bkg against each other in the side band regions """

    mids = [(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]
    x_err = (bins[1]-bins[0])/2.
    fig, ax = plt.subplots(1,1)
    ax.errorbar(mids, data, label='data', 
                color='black', marker='o', linestyle=None, capsize=0, capthick=0)
    ax.errorbar(mids, bkg, label='background', 
                color='red', capsize=0., capthick=0, )

    ax.legend(loc='best')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    fig.set_size_inches(16, 9)
    fig.savefig('plots/'+title+'.png')
    plt.close(fig)


def main():
    """ derives the scale of background using data sidebands """
    # options
    parser = ap.ArgumentParser(description="options for bkg scale derivation")
    parser.add_argument("-c", "--config", default=None, type=str,
                        help="config file for bkg scale derivation")
    parser.add_argument("--lumi-scale", dest='lumi_scale', default=41.5, 
                        type=float, help="luminosity scale for bkg")

    args = parser.parse_args()

    dat_df, bkg_df = load_cfg_to_df(args.config)

    dat_side_bands, bins = get_side_bands_hist(dat_df, 'CMS_hgg_mass', _kIsData=True)
    bkg_side_bands, bins = get_side_bands_hist(bkg_df, 'CMS_hgg_mass', lumi_scale=args.lumi_scale)

    plot(bins, 
        dat_side_bands, bkg_side_bands, 
        "Diphoton Invariant Mass [GeV]", "Events / 0.05 GeV ",
        "mass_sidebands")

    dat_dipho_mva, bins = get_side_bands_hist(dat_df, 'diphoton_mva', hrange=[0,1], num_bins=50, _kIsData=True, lumi_scale=1.)
    bkg_dipho_mva, bins = get_side_bands_hist(bkg_df, 'diphoton_mva', hrange=[0,1], num_bins=50, lumi_scale=args.lumi_scale)

    plot(bins, 
        dat_dipho_mva, bkg_dipho_mva, 
        "Diphoton MVA Score", "Events / 0.02",
        "diphoton_mva_sidebands")

    events_data = np.sum(dat_side_bands)
    events_bkg = np.sum(bkg_side_bands)

    print(f'there are {events_data} in the side bands in data')
    print(f'there are {events_bkg} in the side bands in bkg')
    print(f'the scale for bkg is: {events_data/events_bkg}')


if __name__ == "__main__":
    main()
