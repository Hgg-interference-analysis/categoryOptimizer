import numpy as np
import matplotlib.pyplot as plt
import uproot3 as up
import pandas as pd

from python.helpers.helper_plot import color_gradient

def main():

    LUMI_2017 = 41.5
    MC_TO_DD = 2.7312491505102288

    bkg_files_dd = [
        # ("tagsDumper/trees/Data_13TeV_UntaggedTag_0", "./rootFiles/output_data_UL2017.root", "UL17 Data Sidebands")]
        # ("pp", "./rootFiles/data/ul16_preVFP/ul16_preVFP_prompt-prompt_newBDT.root", "UL16 PreVFP prompt-rompt"),
        # ("pp", "./rootFiles/data/ul16_postVFP/ul16_postVFP_prompt-prompt_newBDT.root", "UL16 PostVFP prompt-rompt"),
        # ("DataDriven_QCD", "./rootFiles/data/ul16_preVFP/ul16_preVFP_prompt-fake_fake-fake_newBDT.root", "UL16 PreVFP prompt-rompt"),
        # ("DataDriven_QCD", "./rootFiles/data/ul16_postVFP/ul16_postVFP_prompt-fake_fake-fake_newBDT.root", "UL16 PostVFP prompt-rompt"),
    ]
    bkg_files_mc = [
        # ("pp","./rootFiles/ul18_dataDriven_newBDT_prompt-prompt_bkg.root", "prompt-prompt"),
        # ("DataDriven_QCD","./rootFiles/ul18_dataDriven_newBDT_prompt-fake_fake-fake_bkg.root", "prompt-fake + fake-fake"),
        #("pp","./rootFiles/data_driven_prompt-prompt_ul17.root", "prompt-prompt"),
        #("DataDriven_QCD","./rootFiles/data_driven_prompt-fake_fake-fake_ul17.root", "prompt-fake + fake-fake"),
        # ("pp", "./rootFiles/data/ul16_preVFP/ul16_preVFP_prompt-prompt_newBDT.root", "UL16 PreVFP prompt-prompt"),
        # ("DataDriven_QCD", "./rootFiles/data/ul16_preVFP/ul16_preVFP_prompt-fake_fake-fake_newBDT.root", "UL16 PreVFP prompt-fake + fake-fake"),
        ("pp", "./rootFiles/data/ul16_postVFP/ul16_postVFP_prompt-prompt_newBDT.root", "UL16 PostVFP prompt-prompt"),
        ("DataDriven_QCD", "./rootFiles/data/ul16_postVFP/ul16_postVFP_prompt-fake_fake-fake_newBDT.root", "UL16 PostVFP prompt-fake + fake-fake"),
    ]

    sig_files_mc = [
        #("tagsDumper/trees/ggh_125_13TeV_UntaggedTag_0","./rootFiles/ul17_ggh_oldBDT.root","ul17 ggh_M125"),
        #("tagsDumper/trees/ggh_125_13TeV_UntaggedTag_0","./rootFiles/ul17_newDMVA_ggh.root","ul17 ggh_M125"),
        #("tagsDumper/trees/wh_125_13TeV_UntaggedTag_0","./rootFiles/ul17_newDMVA_wh.root","ul17 vh_M125"),
        #("tagsDumper/trees/vbf_125_13TeV_UntaggedTag_0","./rootFiles/ul17_newDMVA_vbf.root","ul17 vbf_M125"),
        #("tagsDumper/trees/thq_125_13TeV_UntaggedTag_0","./rootFiles/ul17_newDMVA_thq.root","ul17 thq_M125"),
        #("tagsDumper/trees/thw_125_13TeV_UntaggedTag_0","./rootFiles/ul17_newDMVA_thw.root","ul17 thw_M125"),
        #("tagsDumper/trees/tth_125_13TeV_UntaggedTag_0","./rootFiles/ul17_newDMVA_tth.root","ul17 ttH_M125"),
        #("tagsDumper/trees/ggh_125_13TeV_UntaggedTag_0","./rootFiles/ggH_125_ul18.root","ul18 ggh_M125"),
        #("tagsDumper/trees/wh_125_13TeV_UntaggedTag_0","./rootFiles/wh_125_ul18.root","ul18 vh_M125"),
        #("tagsDumper/trees/vbf_125_13TeV_UntaggedTag_0","./rootFiles/vbf_125_ul18.root","ul18 vbf_M125"),
        #("tagsDumper/trees/tth_125_13TeV_UntaggedTag_0","./rootFiles/ttH_125_ul18.root","ul18 ttH_M125"),
        #("tagsDumper/trees/ggh_125_13TeV_UntaggedTag_0","./rootFiles/ul18_newDMVA_ggh.root","ul18 ggh_M125"),
        #("tagsDumper/trees/wh_125_13TeV_UntaggedTag_0","./rootFiles/ul18_newDMVA_wh.root","ul18 vh_M125"),
        #("tagsDumper/trees/vbf_125_13TeV_UntaggedTag_0","./rootFiles/ul18_newDMVA_vbf.root","ul18 vbf_M125"),
        #("tagsDumper/trees/tth_125_13TeV_UntaggedTag_0","./rootFiles/ul18_newDMVA_tth.root","ul18 ttH_M125"),
        # ("Sig125", "./rootFiles/data/ul16_preVFP/ul16_preVFP_signal_newBDT.root", "UL16 PreVFP Signal"),
        ("Sig125", "./rootFiles/data/ul16_postVFP/ul16_postVFP_signal_newBDT.root", "UL16 PostVFP Signal"),
        # ("Sig125", "./rootFiles/ul18_signal_newBDT.root", "UL18 Signal"),
    ]


    keep = ['CMS_hgg_mass','diphoton_mva','weight']

    df_bkg_dd = []
    df_bkg_mc = []
    df_sig_mc = []

    for f in bkg_files_dd:
        df = up.open(f[1])[f[0]].pandas.df(keep)
        df_bkg_dd.append(df)

    for f in bkg_files_mc:
        df_bkg_mc.append(up.open(f[1])[f[0]].pandas.df(keep))

    for f in sig_files_mc:
        df_sig_mc.append(up.open(f[1])[f[0]].pandas.df(keep))

    hists_bkg_dd = []
    hists_bkg_mc = []
    hists_sig_mc = []

    bins = []
    n_bins = 100
    for df in df_bkg_dd:
        h, bins = np.histogram(df['diphoton_mva'].values, bins=n_bins, range=(0,1) )
        hists_bkg_dd.append(h)

    for i,df in enumerate(df_bkg_mc):
        scale = 1.
        h, bins = np.histogram(df['diphoton_mva'].values, weights = scale*df['weight'].values, bins=n_bins, range=(0,1) )
        hists_bkg_mc.append(h)

    for df in df_sig_mc:
        scale = LUMI_2017
        h, bins = np.histogram(df['diphoton_mva'].values, weights = scale*df['weight'].values, bins=n_bins, range=(0,1) )
        hists_sig_mc.append(h)

    colors = color_gradient(len(bkg_files_dd)+len(bkg_files_mc)+len(sig_files_mc))
    mids = [(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]

    fig, ax = plt.subplots(1,1)

    
    # for i in range(len(hists_bkg_dd)):
    #     stack = hists_bkg_dd[0]
    #     for j in range(i+1, len(hists_bkg_dd)):
    #         stack = np.add(stack, hists_bkg_dd[j])
    #     stack_full = [stack[0]]+list(stack)+[stack[-1]]
    #     ax.plot(mids, stack, 'k.', label=f"{bkg_files_dd[i][2]}")

    
    for i in range(len(hists_bkg_mc)):
        stack = hists_bkg_mc[i]
        for j in range(i+1, len(hists_bkg_mc)):
            stack = np.add(stack, hists_bkg_mc[j])
        stack_full = [stack[0]]+list(stack)+[stack[-1]]
        ax.fill_between([bins[0]]+mids+[bins[-1]], stack_full, color=colors[len(hists_bkg_dd)+i], step="mid", label=f"{bkg_files_mc[i][2]}",)
        ax.errorbar(mids, stack,
                    xerr=(bins[1]-bins[0])/2,
                    yerr=0,
                    drawstyle="steps-mid", color="black", linewidth=0.5)
    

    for i in range(len(hists_sig_mc)):
        stack = hists_sig_mc[i]
        for j in range(i+1, len(hists_sig_mc)):
            stack = np.add(stack, hists_sig_mc[j])
        stack_full = [stack[0]]+list(stack)+[stack[-1]]
        ax.fill_between([bins[0]]+mids+[bins[-1]], stack_full, color=colors[len(hists_bkg_dd)+len(hists_bkg_mc)+i], step="mid", label=f"1000*{sig_files_mc[i][2]}",)
        ax.errorbar(mids, stack,
                    xerr=(bins[1]-bins[0])/2,
                    yerr=0,
                    drawstyle="steps-mid", color="black", linewidth=0.5)

    ax.set_yscale('log')
    ax.set_ybound(10**(-2), 10**7)
    ax.legend()
    ax.set_ylabel('Event/0.01')
    ax.set_xlabel('Diphoton MVA Score')

    plt.show()

    integral_dd = sum(hists_bkg_dd[i][j] for i in range(len(hists_bkg_dd)) for j in range(len(hists_bkg_dd[i])))
    integral_mc = sum(hists_bkg_mc[i][j] for i in range(len(hists_bkg_mc)) for j in range(len(hists_bkg_mc[i])))
    integral_signal = sum(hists_sig_mc[i][j] for i in range(len(hists_sig_mc)) for j in range(len(hists_sig_mc[i])))

    print(integral_dd, integral_mc, integral_dd/integral_mc)
    print(integral_signal)



if __name__=='__main__':
    main()
