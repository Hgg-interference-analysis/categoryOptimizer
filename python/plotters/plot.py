import numpy as np
import pandas as pd

import python.helpers.helper_plot as plot


def xcheck_plot(df, output_tag):
    """
    plot a bunch of cross checks:
    R9 and pt distributions for signal and background
    """

    do_not_plot = ['run','candidate_id', 'weight', 'nvtx', 'processIndex']
    for col in df.columns:
        if col in do_not_plot:
            continue
        plot.plot_df_var(col, df[col].values, df['weight'].values, output_tag)
    
    return


def stack_plot(sig, bkg, output_tag, bounds=[]):
    """ makes a stacked histogram plot """

    lumi_scale = 41.5 # 2017 luminosity
    bkg_scale = 2.76  # side-band estimate of num bkg events

    sig_dfs, sig_titles = sig
    bkg_dfs, bkg_titles = bkg

    labels = list(sig_titles) + list(bkg_titles)

    # collect signal and background
    hists_sig, err_sig, bins = plot.collect_hists(sig_dfs, lumi_scale=lumi_scale*100)
    hists_bkg, err_bkg, bins = plot.collect_hists(bkg_dfs, lumi_scale=lumi_scale, bkg_scale=bkg_scale)

    # create stack plots
    stack_hist = np.zeros(len(bins)-1)
    hists_bkg = hists_bkg[::-1]
    for i,hist in enumerate(hists_bkg):
        stack_hist = np.add(stack_hist, hist)
        hists_bkg[i] = stack_hist

    # setup the plot
    hists_bkg = hists_bkg[::-1]
    hists = hists_sig + hists_bkg
    bin_errors = err_sig + err_bkg
    
    plot.plot_stacked_hists(hists, bins, bin_errors, labels, output_tag, bounds)