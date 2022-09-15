import numpy as np
import matplotlib.pyplot as plt
import uproot3 as up
import pandas as pd

from python.helpers.helper_plot import color_gradient

def main():

    LUMI_2017 = 41.5

    bkg_files_dd = [
        ("pp", "./rootFiles/data_driven_prompt-prompt_ul17.root", "data driven prompt-prompt"),
        ("DataDriven_QCD","./rootFiles/data_driven_prompt-fake_fake-fake_ul17.root", "data driven prompt-fake + fake-fake")
    ]

    bkg_files_mc = [
        ("tagsDumper/trees/prompt_prompt_13TeV_UntaggedTag_0","./rootFiles/gg_bkg_ul17.root", "MC prompt-prompt"),
        ("tagsDumper/trees/prompt_fake_13TeV_UntaggedTag_0","./rootFiles/gjet_bkg_ul17.root", "MC prompt-fake")
    ]

    sig_files_mc = [
        ("tagsDumper/trees/tth_125_13TeV_UntaggedTag_0","./rootFiles/tth_125_ul17.root", "ttH"),
        ("tagsDumper/trees/vbf_125_13TeV_UntaggedTag_0","./rootFiles/vbf_125_ul17.root", "VBF"),
        ("tagsDumper/trees/thw_125_13TeV_UntaggedTag_0","./rootFiles/thw_125_ul17.root", "tHW"),
        ("tagsDumper/trees/thq_125_13TeV_UntaggedTag_0","./rootFiles/thq_125_ul17.root", "tHq"),
        ("tagsDumper/trees/vh_125_13TeV_UntaggedTag_0","./rootFiles/vh_125_ul17.root", "VH"),
        ("tagsDumper/trees/ggh_125_13TeV_UntaggedTag_0","./rootFiles/ggh_125_ul17.root", "ggH")
    ]

    keep = ['CMS_hgg_mass','diphoMVA','weight']

    df_bkg_dd = []
    df_bkg_mc = []
    df_sig_mc = []

    for f in bkg_files_dd:
        df_bkg_dd.append(up.open(f[1])[f[0]].pandas.df(keep))

    for f in bkg_files_mc:
        df_bkg_mc.append(up.open(f[1])[f[0]].pandas.df(keep))

    for f in sig_files_mc:
        df_sig_mc.append(up.open(f[1])[f[0]].pandas.df(keep))

    hists_bkg_dd = []
    hists_bkg_mc = []
    hists_sig_mc = []

    bins = []
    for df in df_bkg_dd:
        h, bins = np.histogram(df['diphoMVA'].values, weights = df['weight'].values, bins=200, range=(0,1) )
        hists_bkg_dd.append(h)
    
    for df in df_bkg_mc:
        h, bins = np.histogram(df['diphoMVA'].values, weights = LUMI_2017*df['weight'].values, bins=200, range=(0,1) )
        hists_bkg_mc.append(h)

    for df in df_sig_mc:
        h, bins = np.histogram(df['diphoMVA'].values, weights = df['weight'].values, bins=200, range=(0,1) )
        hists_sig_mc.append(h)

    colors = color_gradient(len(bkg_files_dd)+len(bkg_files_mc)+len(sig_files_mc))

    mids = [(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]
    
    fig, ax = plt.subplots(1,1)

    for i in range(len(hists_bkg_dd)):
        stack = hists_bkg_dd[0]
        for j in range(i+1, len(hists_bkg_dd)):
            stack = np.add(stack, hists_bkg_dd[j])
        stack_full = [stack[0]]+list(stack)+[stack[-1]]
        ax.fill_between([bins[0]]+mids+[bins[-1]], stack_full, color=colors[i])
        ax.errorbar(mids, stack, 
                    xerr=(bins[1]-bins[0])/2,
                    yerr=np.sqrt(stack),
                    label=f"{bkg_files_dd[i][2]}",
                    drawstyle="steps-mid", color=colors[i])

    for i in range(len(hists_bkg_mc)):
        stack = hists_bkg_mc[0]
        for j in range(i+1, len(hists_bkg_mc)):
            stack = np.add(stack, hists_bkg_mc[j])
        stack_full = [stack[0]]+list(stack)+[stack[-1]]
        ax.fill_between([bins[0]]+mids+[bins[-1]], stack_full, color=colors[len(hists_bkg_dd)+i])
        ax.errorbar(mids, stack, 
                    xerr=(bins[1]-bins[0])/2,
                    yerr=np.sqrt(stack),
                    label=f"{bkg_files_dd[i][2]}",
                    drawstyle="steps-mid", color=colors[len(hists_bkg_dd)+i]) 

    for i in range(len(hists_sig_mc)):
        stack = hists_sig_mc[0]
        for j in range(i+1, len(hists_sig_mc)):
            stack = np.add(stack, hists_sig_mc[j])
        stack_full = [stack[0]]+list(stack)+[stack[-1]]
        ax.fill_between([bins[0]]+mids+[bins[-1]], stack_full, color=colors[len(hists_bkg_dd)+len(hists_bkg_mc)+i])
        ax.errorbar(mids, stack, 
                    xerr=(bins[1]-bins[0])/2,
                    yerr=np.sqrt(stack),
                    label=f"{bkg_files_dd[i][2]}",
                    drawstyle="steps-mid", color=colors[len(hists_bkg_dd)+len(hists_bkg_mc)+i])

    plt.show()


if __name__=='__main__':
    main()    