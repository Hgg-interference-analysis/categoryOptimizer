import numpy as np
from numpy.core.shape_base import stack
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import iqr


def color_gradient(n):
    """ returns a list of rgb values scaled from blue to green to red """

    colors = []
    for i in range(n):
        blue = (1 - 2*float(i)/float(n-1))
        if blue < 0:
            blue = 0.
        red = float(i)/float(n-1) if blue == 0 else 0.
        green = 0.0
        if i == (n-1)/2:
            green = 1.
        elif i < (n-1)/2:
            green = 2 * float(i)/float(n-1)
        else:
            green = 2 - 2*float(i)/float(n-1)

        colors.append((red, green, blue))

    return colors


def get_bin_uncertainties(bins, values, weights):
    """ calculates the uncertainty of weighted bins """

    ret = []
    for i in range(len(bins)-1):
        val_mask = np.logical_and(bins[i] <= values, values < bins[i+1])
        ret.append(np.sqrt(np.sum(np.power(weights[val_mask], 2))))

    return np.array(ret)


def plot_df_var(var_name, values, weights, output_tag):
    """ plots var_name as a histogram """

    plot_name = f'./plots/{output_tag}_{var_name}.png'
    print(f'[INFO] plotting {var_name} in {plot_name}')
    values = np.array(np.abs(values))
    weights = np.array(weights)
    bin_width = 2*iqr(np.multiply(values, weights)) / \
        np.power(np.sum(np.multiply(values, weights)), 1/3)
    if bin_width < 0.1:
        bin_width = 0.1
    num_bins = int((max(values)-min(values))//bin_width)
    hist, bins = np.histogram(values, bins=num_bins, weights=weights)
    mids = [(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]
    mids_full = mids.copy()
    hist_full = list(hist.copy())
    hist_full.append(hist[-1])
    mids_full.append(bins[-1])
    mids_full.insert(0, bins[0])
    hist_full.insert(0, hist[0])

    x_err = (bins[1]-bins[0])/2
    y_err = get_bin_uncertainties(bins, values, weights)

    fig, axs = plt.subplots(nrows=1, ncols=1)
    fig.subplots_adjust(left=0.12, right=0.88, top=0.95, bottom=0.1)

    axs.fill_between(mids_full, hist_full,
                     y2=0, step='mid',
                     alpha=0.2, color='cornflowerblue')
    axs.errorbar(mids, hist,
                 xerr=x_err, yerr=y_err,
                 label=f'{var_name}',
                 drawstyle='steps-mid', capsize=0.)

    axs.legend(loc='best')
    axs.set_ylabel(f'Events / {2*x_err}')
    axs.set_xlabel(f'{var_name}')

    plt.grid(which='major', axis='both')
    fig.savefig(plot_name)
    plt.close(fig)


def collect_hists(dfs, lumi_scale=1., bkg_scale=1.):
    """ returns histograms from a list of dataframes """

    hists = []
    bin_errors = []
    bins = []
    num_bins = 50

    for i, df in enumerate(dfs):
        hist, bins = np.histogram(df['diphoton_mva'].values, bins=num_bins, range=[
                                  0, 1], weights=df['weight'].values)
        hist = lumi_scale * hist
        bin_errors.append(get_bin_uncertainties(
            bins, df['diphoton_mva'].values, df['weight'].values))
        hists.append(np.copy(hist))

    return hists, bin_errors, bins


def plot_stacked_hists(hists, bins, bin_errors, labels, output_tag, bounds):
    """ creates a stacked histogram plot """

    stack_colors = color_gradient(len(hists))
    fig, axs = plt.subplots(1, 1)
    fig.subplots_adjust(left=0.085, right=0.99, top=0.95, bottom=0.1)

    mids = [(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]
    mids_full = mids.copy()
    mids_full.append(bins[-1])
    mids_full.insert(0, bins[0])

    x_err = (bins[1]-bins[0])/2

    # plot the histograms
    handels = []
    for i, hist in enumerate(hists):
        hist_full = list(hist).copy()
        hist_full.insert(0, hist_full[0])
        hist_full.insert(-1, hist_full[-1])
        next_hist = 0
        if i < len(hists)-1:
            next_hist = list(hists[i+1].copy())
            next_hist.insert(0, next_hist[0])
            next_hist.insert(-1, next_hist[-1])

        #axs.fill_between(mids_full, hist_full,
        #                 y2=next_hist, step='mid',
        #                 alpha=1, color=stack_colors[::-1][i])

        errorbar = axs.errorbar(mids, hist,
                                xerr=x_err, yerr=bin_errors[i],
                                drawstyle='steps-mid', capsize=0., 
                                linewidth=3,
                                color=stack_colors[::-1][i])
        fill = axs.fill(np.NaN, np.NaN, 
                        color=stack_colors[::-1][i], alpha=0.5)
        handels.append((fill[0], errorbar[0]))

    for i,bound in enumerate(bounds):
        axs.plot([bound, bound], [0.3, 10**7], 
                color='black', linestyle=(0, (5, 5)), markersize=0)
        axs.annotate(f"Untagged {len(bounds)-1-i}", (bound, 0.3), xycoords='data',
                    ha='left', va='bottom', rotation=270)

    axs.legend(handels, [f'{x}' for x in labels], loc='upper right')

    axs.set_xlabel("Diphoton MVA Score", ha='right', x=1.)
    axs.set_ylabel("Events / 0.02", ha='right', y=1.)
    plt.ylim(0.3, 10**7)
    axs.set_yscale('log')
    axs.annotate("$\\bf{CMS} \ \\it{Simulation \ Preliminary}$",
                xy=(0,1), xycoords="axes fraction",
                ha='left', va='bottom')
    axs.annotate('41.5 fb$^{-1}$ (13 TeV) 2017',
                xy=(1,1), xycoords="axes fraction",
                ha='right', va='bottom')

    plt.grid(which='major', axis='both')
    fig.set_size_inches(7, 5)
    fig.savefig('plots/stacked_histogram_'+output_tag+'.png')
    plt.close(fig)
