import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import iqr

def get_bin_uncertainties(bins, values, weights):
    """ calculates the uncertainty of weighted bins """
    ret = []
    for i in range(len(bins)-1):
        val_mask = np.logical_and( bins[i] <= values, values < bins[i+1])
        ret.append(np.sqrt(np.sum(np.power(weights[val_mask], 2))))

    return np.array(ret)


def plot_df_var(var_name, values, weights, output_tag):
    """ plots var_name as a histogram """
    plot_name = f'./plots/{output_tag}_{var_name}.png'
    print(f'[INFO] plotting {var_name} in {plot_name}')
    values = np.array(np.abs(values))
    weights = np.array(weights)
    bin_width = 2*iqr(np.multiply(values, weights))/np.power(np.sum(np.multiply(values,weights)),1/3)
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
    
    fig,axs = plt.subplots(nrows=1, ncols=1)
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


def plot(df, output_tag):
    """
    plot a bunch of cross checks:
    R9 and pt distributions for signal and background
    """

    do_not_plot = ['run','candidate_id', 'weight', 'nvtx', 'processIndex']
    for col in df.columns:
        if col in do_not_plot:
            continue
        plot_df_var(col, df[col].values, df['weight'].values, output_tag)
    
    return
