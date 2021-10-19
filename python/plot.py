import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc

def plot_unstacked(files):
    """
    files is composed of two signal and two background samples
    """
    print("plotting unstacked plots")
    
    dat_sig = pd.concat(files[0:2])
    dat_bkg = pd.concat(files[2::])

    sig_r9, bins_r9 = np.histogram(dat_sig['lead_R9'].values, bins=100, range=(0.5,1.1))
    sig_subleadr9, bins_r9 = np.histogram(dat_sig['sublead_R9'].values, bins=100, range=(0.5,1.1))
    sig_r9 = np.add(sig_r9,sig_subleadr9)
    sig_r9 = sig_r9/np.sum(sig_r9)
    bkg_r9, bins_r9 = np.histogram(dat_bkg['lead_R9'].values, bins=100, range=(0.5,1.1))
    bkg_subr9, bins_r9 = np.histogram(dat_bkg['sublead_R9'].values, bins=100, range=(0.5,1.1))
    bkg_r9 = np.add(bkg_r9,bkg_subr9)
    bkg_r9 = bkg_r9/np.sum(bkg_r9)

    sig_pt, bins_pt = np.histogram(dat_sig['leadPt'].values, bins=400, range=(0,200))
    sig_subpt, bins_pt = np.histogram(dat_sig['subleadPt'].values, bins=400, range=(0,200))
    sig_pt = np.add(sig_pt, sig_subpt)
    sig_pt = sig_pt/np.sum(sig_pt)
    bkg_pt, bins_pt = np.histogram(dat_bkg['leadPt'].values, bins=400, range=(0,200))
    bkg_subpt, bins_pt = np.histogram(dat_bkg['subleadPt'].values, bins=400, range=(0,200))
    bkg_pt = np.add(bkg_pt,bkg_subpt)
    bkg_pt = bkg_pt/np.sum(bkg_pt)
        
    fig, axs = plt.subplots(nrows=1, ncols=1)
    axs.set_xlim(0.5,1.5)
    axs.set_ylim(0, 1.2*max(sig_r9.max(),bkg_r9.max()))
    axs.hist(bins_r9[:-1]+0.5*(bins_r9[1:]-bins_r9[:-1]),weights=sig_r9, bins=bins_r9, label='Signal', alpha=0.5)
    axs.hist(bins_r9[:-1]+0.5*(bins_r9[1:]-bins_r9[:-1]),weights=bkg_r9, bins=bins_r9, label='Background', alpha=0.5)
    axs.set_xlabel("R$_{9}$")
    axs.set_ylabel("Fraction")
    axs.legend(loc='best')
    fig.savefig('/afs/cern.ch/work/n/nschroed/ss_pyfit/CMSSW_10_2_14/src/diphoton-mva-optimizer/plots/dists_r9.png')
    plt.close(fig)

    fig, axs = plt.subplots(nrows=1, ncols=1)
    axs.set_xlim(0,200)
    axs.set_ylim(0, 1.2*max(sig_pt.max(),bkg_pt.max()))
    axs.hist(bins_pt[:-1]+0.5*(bins_pt[1:]-bins_pt[:-1]),weights=sig_pt, bins=bins_pt, label='Signal', alpha=0.5)
    axs.hist(bins_pt[:-1]+0.5*(bins_pt[1:]-bins_pt[:-1]),weights=bkg_pt, bins=bins_pt, label='Background', alpha=0.5)
    axs.set_xlabel("P$_{T}$ [GeV]")
    axs.set_ylabel("Fraction / 0.5 GeV")
    axs.legend(loc='best')
    fig.savefig('/afs/cern.ch/work/n/nschroed/ss_pyfit/CMSSW_10_2_14/src/diphoton-mva-optimizer/plots/dists_pt.png')
    plt.close(fig)

def plot_stacked(files):
    """
    Stacked histograms for diphotonMVA, leading and subleading idMVA,
    and sigmaMDecorr for background and combined signal
    DiphotonMVA
    decorrSigmaM
    leadIDMVA
    subleadIDMVA
    """
    print("plotting stacked plots")

    #bkg0 = diphoJetsBox
    #bkg1 = gjets
    dat_sig = pd.concat(files[0:2])
    dat_bkg0 = files[2]
    dat_bkg1 = files[3]

    lumi = 41.5

    mask_bkg0_promptPrompt = (dat_bkg0['gen_lead_pt'] != 0).values & (dat_bkg0['gen_sublead_pt'] != 0).values
    mask_bkg0_promptFake = ((dat_bkg0['gen_lead_pt'] != 0).values & (dat_bkg0['gen_sublead_pt'] == 0).values) | ((dat_bkg0['gen_lead_pt'] == 0).values & (dat_bkg0['gen_sublead_pt'] != 0).values) 
    mask_bkg0_fakeFake = (dat_bkg0['gen_lead_pt'] == 0).values & (dat_bkg0['gen_sublead_pt'] == 0).values

    mask_bkg1_promptPrompt = (dat_bkg1['gen_lead_pt'] != 0).values & (dat_bkg1['gen_sublead_pt'] != 0).values
    mask_bkg1_promptFake = ((dat_bkg1['gen_lead_pt'] != 0).values & (dat_bkg1['gen_sublead_pt'] == 0).values) | ((dat_bkg1['gen_lead_pt'] == 0).values & (dat_bkg1['gen_sublead_pt'] != 0).values) 
    mask_bkg1_fakeFake = (dat_bkg1['gen_lead_pt'] == 0).values & (dat_bkg1['gen_sublead_pt'] == 0).values

    mask_bkg0_leadFake = dat_bkg0['gen_lead_pt'] == 0
    mask_bkg1_leadFake = dat_bkg1['gen_lead_pt'] == 0

    mask_bkg0_subleadFake = dat_bkg0['gen_sublead_pt'] == 0
    mask_bkg1_subleadFake = dat_bkg1['gen_sublead_pt'] == 0


    dat_bkg0 = dat_bkg0[mask_bkg0_leadFake]
    dat_bkg1 = dat_bkg1[mask_bkg1_leadFake]

    print("diphoton jets box breakdown:")
    print("prompt-prompt:", np.sum(mask_bkg0_promptPrompt))
    print("prompt-fake:", np.sum(mask_bkg0_promptFake))
    print("fake-fake:", np.sum(mask_bkg0_fakeFake))

    print("gamma plus jets breakdown:")
    print("prompt-prompt:", np.sum(mask_bkg1_promptPrompt))
    print("prompt-fake:", np.sum(mask_bkg1_promptFake))
    print("fake-fake:", np.sum(mask_bkg1_fakeFake))

    sig_diphoMVA, bins_dipho = np.histogram(dat_sig['DiphotonMVA'].values, bins=50, range=(0,1), weights=dat_sig['weight'].values)
    sig_decorr, bins_decorr = np.histogram(dat_sig['decorrSigmaM'].values, bins=100, range=(0,0.05), weights=dat_sig['weight'].values)
    sig_lead_idmva, bins_idmva = np.histogram(dat_sig['leadIDMVA'].values, bins=200, range=(-1,1), weights=dat_sig['weight'].values)
    sig_sublead_idmva, bins_idmva = np.histogram(dat_sig['subleadIDMVA'].values, bins=200, range=(-1,1), weights=dat_sig['weight'].values)
    sig_min_idmva, bins_idmva = np.histogram(np.minimum(dat_sig['leadIDMVA'].values, dat_sig['subleadIDMVA'].values), bins=200, range=(-1,1), weights=dat_sig['weight'].values)

    sig_diphoMVA = sig_diphoMVA*lumi
    sig_decorr = sig_decorr*lumi
    sig_lead_idmva = sig_lead_idmva*lumi
    sig_sublead_idmva = sig_sublead_idmva*lumi
    sig_min_idmva = sig_min_idmva*lumi

    bkg0_diphoMVA, bins_dipho = np.histogram(dat_bkg0['DiphotonMVA'].values, bins=50, range=(0,1), weights=dat_bkg0['weight'].values)
    bkg0_decorr, bins_decorr = np.histogram(dat_bkg0['decorrSigmaM'].values, bins=100, range=(0,0.05), weights=dat_bkg0['weight'].values)
    bkg0_lead_idmva, bins_idmva = np.histogram(dat_bkg0['leadIDMVA'].values, bins=200, range=(-1,1), weights=dat_bkg0['weight'].values)
    bkg0_sublead_idmva, bins_idmva = np.histogram(dat_bkg0['subleadIDMVA'].values, bins=200, range=(-1,1), weights=dat_bkg0['weight'].values)
    bkg0_min_idmva, bins_idmva = np.histogram(np.minimum(dat_bkg0['leadIDMVA'].values, dat_bkg0['subleadIDMVA'].values), bins=200, range=(-1,1), weights=dat_bkg0['weight'].values)

    bkg0_diphoMVA = bkg0_diphoMVA*lumi
    bkg0_decorr = bkg0_decorr*lumi
    bkg0_lead_idmva = bkg0_lead_idmva*lumi
    bkg0_sublead_idmva = bkg0_sublead_idmva*lumi
    bkg0_min_idmva = bkg0_min_idmva*lumi

    bkg1_diphoMVA, bins_dipho = np.histogram(dat_bkg1['DiphotonMVA'].values, bins=50, range=(0,1), weights=dat_bkg1['weight'].values)
    bkg1_decorr, bins_decorr = np.histogram(dat_bkg1['decorrSigmaM'].values, bins=100, range=(0,0.05), weights=dat_bkg1['weight'].values)
    bkg1_lead_idmva, bins_idmva = np.histogram(dat_bkg1['leadIDMVA'].values, bins=200, range=(-1,1), weights=dat_bkg1['weight'].values)
    bkg1_sublead_idmva, bins_idmva = np.histogram(dat_bkg1['subleadIDMVA'].values, bins=200, range=(-1,1), weights=dat_bkg1['weight'].values)
    bkg1_min_idmva, bins_idmva = np.histogram(np.minimum(dat_bkg1['leadIDMVA'].values, dat_bkg1['subleadIDMVA'].values), bins=200, range=(-1,1), weights=dat_bkg1['weight'].values)

    bkg1_diphoMVA = bkg1_diphoMVA*lumi
    bkg1_decorr = bkg1_decorr*lumi
    bkg1_lead_idmva = bkg1_lead_idmva*lumi
    bkg1_sublead_idmva = bkg1_sublead_idmva*lumi
    bkg1_min_idmva = bkg1_min_idmva*lumi

    bkg_tot_diphoMVA = np.add(bkg1_diphoMVA,bkg0_diphoMVA)
    bkg_tot_decorr = np.add(bkg1_decorr,bkg0_decorr)
    bkg_tot_lead_idmva = np.add(bkg1_lead_idmva,bkg0_lead_idmva)
    bkg_tot_sublead_idmva = np.add(bkg1_sublead_idmva,bkg0_sublead_idmva)
    bkg_tot_min_idmva = np.add(bkg0_min_idmva,bkg1_min_idmva)

    sigPlusBkg_diphoMVA = np.add(sig_diphoMVA,bkg_tot_diphoMVA)
    sigPlusBkg_decorr = np.add(sig_decorr,bkg_tot_decorr)
    sigPlusBkg_lead_idmva = np.add(sig_lead_idmva,bkg_tot_lead_idmva)
    sigPlusBkg_sublead_idmva = np.add(sig_sublead_idmva,bkg_tot_sublead_idmva)
    sigPlusBkg_min_idmva = np.add(sig_min_idmva,bkg_tot_min_idmva)

    """
    #diphoMVA 

    #stacked background and signal
    fig, axs = plt.subplots(nrows=1, ncols=1)
    axs.set_xlim(0,1)
    axs.set_ylim(0.001, 1.2*max(sigPlusBkg_diphoMVA.max(),bkg_tot_diphoMVA.max()))
    axs.hist(bins_dipho[:-1]+0.5*(bins_dipho[1:]-bins_dipho[:-1]),weights=sigPlusBkg_diphoMVA, bins=bins_dipho, label='Signal')
    axs.hist(bins_dipho[:-1]+0.5*(bins_dipho[1:]-bins_dipho[:-1]),weights=bkg_tot_diphoMVA, bins=bins_dipho, label='GJets')
    axs.hist(bins_dipho[:-1]+0.5*(bins_dipho[1:]-bins_dipho[:-1]),weights=bkg0_diphoMVA, bins=bins_dipho, label='DiphotonJetsBox')
    axs.set_xlabel("Diphoton MVA")
    axs.set_ylabel("Events")
    axs.legend(loc='best')
    fig.savefig('/afs/cern.ch/work/n/nschroed/ss_pyfit/CMSSW_10_2_14/src/diphoton-mva-optimizer/plots/stacked_diphoMVA.png')
    axs.set_ylim(10,120*max(bkg0_diphoMVA.max(),bkg1_diphoMVA.max()))
    plt.yscale("log")
    fig.savefig('/afs/cern.ch/work/n/nschroed/ss_pyfit/CMSSW_10_2_14/src/diphoton-mva-optimizer/plots/stacked_diphoMVA_log.png')
    plt.close(fig)
    print("stacked diphoMVA")

    #unstacked
    fig, axs = plt.subplots(nrows=1,ncols=1)
    axs.set_xlim(0,1)
    axs.set_ylim(0.001,1.2*max(bkg0_diphoMVA.max(),bkg1_diphoMVA.max()))
    mids = bins_dipho[:-1]+0.5*(bins_dipho[1:]-bins_dipho[:-1])
    xloerr = mids - bins_dipho[:len(bins_dipho)-1]
    xhierr = bins_dipho[1:] - mids
    x_error = np.array([xloerr,xhierr])
    axs.errorbar(mids, bkg0_diphoMVA, yerr=np.sqrt(bkg0_diphoMVA), xerr=x_error, label="DiphotonJetsBox", color='dodgerblue', marker=",", mew=0, capsize=0, elinewidth=1,linestyle='')
    axs.errorbar(mids, bkg1_diphoMVA, yerr=np.sqrt(bkg1_diphoMVA), xerr=x_error, label="GJets", color='darkblue', marker=",", mew=0, capsize=0, elinewidth=1,linestyle='')
    axs.errorbar(mids, sig_diphoMVA*500, yerr=np.sqrt(sig_diphoMVA*500), xerr=x_error, label="Signal*500", color='darkred', marker=",", mew=0, capsize=0, elinewidth=1,linestyle='')
    axs.set_xlabel("Diphoton MVA")
    axs.set_ylabel("Events")
    axs.legend(loc='best')
    fig.savefig('/afs/cern.ch/work/n/nschroed/ss_pyfit/CMSSW_10_2_14/src/diphoton-mva-optimizer/plots/unstacked_diphoMVA.png')
    axs.set_ylim(10,120*max(bkg0_diphoMVA.max(),bkg1_diphoMVA.max()))
    plt.yscale("log")
    fig.savefig('/afs/cern.ch/work/n/nschroed/ss_pyfit/CMSSW_10_2_14/src/diphoton-mva-optimizer/plots/unstacked_diphoMVA_log.png')
    plt.close(fig)
    print('unstacked diphoMVA')


    fig, axs = plt.subplots(nrows=1, ncols=1)
    axs.set_xlim(0,0.05)
    axs.set_ylim(0, 1.2*max(sigPlusBkg_decorr.max(),bkg_tot_decorr.max()))
    axs.hist(bins_decorr[:-1]+0.5*(bins_decorr[1:]-bins_decorr[:-1]),weights=sigPlusBkg_decorr, bins=bins_decorr, label='Signal')
    axs.hist(bins_decorr[:-1]+0.5*(bins_decorr[1:]-bins_decorr[:-1]),weights=bkg_tot_decorr, bins=bins_decorr, label='GJets')
    axs.hist(bins_decorr[:-1]+0.5*(bins_decorr[1:]-bins_decorr[:-1]),weights=bkg0_decorr, bins=bins_decorr, label='DiphotonJetsBox')
    axs.set_xlabel("Decorr SigmaM [GeV]")
    axs.set_ylabel("Events")
    axs.legend(loc='best')
    fig.savefig('/afs/cern.ch/work/n/nschroed/ss_pyfit/CMSSW_10_2_14/src/diphoton-mva-optimizer/plots/stacked_decorr.png')
    axs.set_ylim(10,120*max(bkg0_decorr.max(),bkg1_decorr.max()))
    plt.yscale("log")
    fig.savefig('/afs/cern.ch/work/n/nschroed/ss_pyfit/CMSSW_10_2_14/src/diphoton-mva-optimizer/plots/stacked_decorr_log.png')
    plt.close(fig)
    print("decorr")

    fig, axs = plt.subplots(nrows=1,ncols=1)
    axs.set_xlim(0,0.05)
    axs.set_ylim(0.001,1.2*max(bkg0_decorr.max(),bkg1_decorr.max()))
    mids = bins_decorr[:-1]+0.5*(bins_decorr[1:]-bins_decorr[:-1])
    xloerr = mids - bins_decorr[:len(bins_decorr)-1]
    xhierr = bins_decorr[1:] - mids
    x_error = np.array([xloerr,xhierr])
    axs.errorbar(mids, bkg0_decorr, yerr=np.sqrt(bkg0_decorr), xerr=x_error, label="DiphotonJetsBox", color='dodgerblue', marker=",", mew=0, capsize=0, elinewidth=1,linestyle='')
    axs.errorbar(mids, bkg1_decorr, yerr=np.sqrt(bkg1_decorr), xerr=x_error, label="GJets", color='darkblue', marker=",", mew=0, capsize=0, elinewidth=1,linestyle='')
    axs.errorbar(mids, sig_decorr*500, yerr=np.sqrt(sig_decorr*500), xerr=x_error, label="Signal*500", color='darkred', marker=",", mew=0, capsize=0, elinewidth=1,linestyle='')
    axs.set_xlabel("Decorr SigmaM [GeV]")
    axs.set_ylabel("Events")
    axs.legend(loc='best')
    fig.savefig('/afs/cern.ch/work/n/nschroed/ss_pyfit/CMSSW_10_2_14/src/diphoton-mva-optimizer/plots/unstacked_decorr.png')
    axs.set_ylim(10,120*max(bkg0_decorr.max(),bkg1_decorr.max()))
    plt.yscale("log")
    fig.savefig('/afs/cern.ch/work/n/nschroed/ss_pyfit/CMSSW_10_2_14/src/diphoton-mva-optimizer/plots/unstacked_decorr_log.png')
    plt.close(fig)
    print('unstacked decorr')

    fig, axs = plt.subplots(nrows=1, ncols=1)
    axs.set_xlim(-1,1)
    axs.set_ylim(0, 1.2*max(sigPlusBkg_lead_idmva.max(),bkg_tot_lead_idmva.max()))
    axs.hist(bins_idmva[:-1]+0.5*(bins_idmva[1:]-bins_idmva[:-1]),weights=sigPlusBkg_lead_idmva, bins=bins_idmva, label='Signal')
    axs.hist(bins_idmva[:-1]+0.5*(bins_idmva[1:]-bins_idmva[:-1]),weights=bkg_tot_lead_idmva, bins=bins_idmva, label='GJets')
    axs.hist(bins_idmva[:-1]+0.5*(bins_idmva[1:]-bins_idmva[:-1]),weights=bkg0_lead_idmva, bins=bins_idmva, label='DiphotonJetsBox')
    axs.set_xlabel("Leading Photon IDMVA")
    axs.set_ylabel("Events")
    axs.legend(loc='best')
    fig.savefig('/afs/cern.ch/work/n/nschroed/ss_pyfit/CMSSW_10_2_14/src/diphoton-mva-optimizer/plots/dists_lead_idmva.png')
    axs.set_ylim(10,120*max(bkg0_lead_idmva.max(),bkg1_lead_idmva.max()))
    plt.yscale("log")
    fig.savefig('/afs/cern.ch/work/n/nschroed/ss_pyfit/CMSSW_10_2_14/src/diphoton-mva-optimizer/plots/stacked_lead_idmva_log.png')
    plt.close(fig)
    print("lead_idmva")

    """
    fig, axs = plt.subplots(nrows=1,ncols=1)
    axs.set_xlim(-1,1)
    axs.set_ylim(0.001,1.2*max(bkg0_lead_idmva.max(),bkg1_lead_idmva.max()))
    mids = bins_idmva[:-1]+0.5*(bins_idmva[1:]-bins_idmva[:-1])
    xloerr = mids - bins_idmva[:len(bins_idmva)-1]
    xhierr = bins_idmva[1:] - mids
    x_error = np.array([xloerr,xhierr])
    axs.errorbar(mids, bkg0_lead_idmva, yerr=np.sqrt(bkg0_lead_idmva), xerr=x_error, label="DiphotonJetsBox", color='dodgerblue', marker=",", mew=0, capsize=0, elinewidth=1,linestyle='')
    axs.errorbar(mids, bkg1_lead_idmva, yerr=np.sqrt(bkg1_lead_idmva), xerr=x_error, label="GJets", color='darkblue', marker=",", mew=0, capsize=0, elinewidth=1,linestyle='')
    #axs.errorbar(mids, sig_lead_idmva*500, yerr=np.sqrt(sig_lead_idmva*500), xerr=x_error, label="Signal*500", color='darkred', marker=",", mew=0, capsize=0, elinewidth=1,linestyle='')
    axs.set_xlabel("Lead IDMVA")
    axs.set_ylabel("Events")
    axs.legend(loc='best')
    fig.savefig('/afs/cern.ch/work/n/nschroed/ss_pyfit/CMSSW_10_2_14/src/diphoton-mva-optimizer/plots/unstacked_lead_idmva_leadFake.png')
    axs.set_ylim(10,120*max(bkg0_lead_idmva.max(),bkg1_lead_idmva.max()))
    plt.yscale("log")
    fig.savefig('/afs/cern.ch/work/n/nschroed/ss_pyfit/CMSSW_10_2_14/src/diphoton-mva-optimizer/plots/unstacked_lead_idmva_log_leadFake.png')
    plt.close(fig)
    print('unstacked lead_idmva')

    """
    fig, axs = plt.subplots(nrows=1, ncols=1)
    axs.set_xlim(-1,1)
    axs.set_ylim(0, 1.2*max(sigPlusBkg_sublead_idmva.max(),bkg_tot_sublead_idmva.max()))
    axs.hist(bins_idmva[:-1]+0.5*(bins_idmva[1:]-bins_idmva[:-1]),weights=sigPlusBkg_sublead_idmva, bins=bins_idmva, label='Signal')
    axs.hist(bins_idmva[:-1]+0.5*(bins_idmva[1:]-bins_idmva[:-1]),weights=bkg_tot_sublead_idmva, bins=bins_idmva, label='GJets')
    axs.hist(bins_idmva[:-1]+0.5*(bins_idmva[1:]-bins_idmva[:-1]),weights=bkg0_sublead_idmva, bins=bins_idmva, label='DiphotonJetsBox')
    axs.set_xlabel("Subleading Photon IDMVA")
    axs.set_ylabel("Events")
    axs.legend(loc='best')
    fig.savefig('/afs/cern.ch/work/n/nschroed/ss_pyfit/CMSSW_10_2_14/src/diphoton-mva-optimizer/plots/dists_sublead_idmva.png')
    axs.set_ylim(10,120*max(bkg0_sublead_idmva.max(),bkg1_sublead_idmva.max()))
    plt.yscale("log")
    fig.savefig('/afs/cern.ch/work/n/nschroed/ss_pyfit/CMSSW_10_2_14/src/diphoton-mva-optimizer/plots/stacked_sublead_idmva_log.png')
    plt.close(fig)
    print("sublead_idmva")

    fig, axs = plt.subplots(nrows=1,ncols=1)
    axs.set_xlim(-1,1)
    axs.set_ylim(0.001,1.2*max(bkg0_sublead_idmva.max(),bkg1_sublead_idmva.max()))
    mids = bins_idmva[:-1]+0.5*(bins_idmva[1:]-bins_idmva[:-1])
    xloerr = mids - bins_idmva[:len(bins_idmva)-1]
    xhierr = bins_idmva[1:] - mids
    x_error = np.array([xloerr,xhierr])
    axs.errorbar(mids, bkg0_sublead_idmva, yerr=np.sqrt(bkg0_sublead_idmva), xerr=x_error, label="DiphotonJetsBox", color='dodgerblue', marker=",", mew=0, capsize=0, elinewidth=1,linestyle='')
    axs.errorbar(mids, bkg1_sublead_idmva, yerr=np.sqrt(bkg1_sublead_idmva), xerr=x_error, label="GJets", color='darkblue', marker=",", mew=0, capsize=0, elinewidth=1,linestyle='')
    #axs.errorbar(mids, sig_sublead_idmva*500, yerr=np.sqrt(sig_sublead_idmva*500), xerr=x_error, label="Signal*500", color='darkred', marker=",", mew=0, capsize=0, elinewidth=1,linestyle='')
    axs.set_xlabel("Subleading IDMVA")
    axs.set_ylabel("Events")
    axs.legend(loc='best')
    fig.savefig('/afs/cern.ch/work/n/nschroed/ss_pyfit/CMSSW_10_2_14/src/diphoton-mva-optimizer/plots/unstacked_sublead_idmva_promptFake.png')
    axs.set_ylim(10,120*max(bkg0_sublead_idmva.max(),bkg1_sublead_idmva.max()))
    plt.yscale("log")
    fig.savefig('/afs/cern.ch/work/n/nschroed/ss_pyfit/CMSSW_10_2_14/src/diphoton-mva-optimizer/plots/unstacked_sublead_idmva_log_promptFake.png')
    plt.close(fig)
    print('unstacked sublead_idmva')

    fig, axs = plt.subplots(nrows=1, ncols=1)
    axs.set_xlim(-1,1)
    axs.set_ylim(0, 1.2*max(sigPlusBkg_min_idmva.max(),bkg_tot_min_idmva.max()))
    axs.hist(bins_idmva[:-1]+0.5*(bins_idmva[1:]-bins_idmva[:-1]),weights=sigPlusBkg_min_idmva, bins=bins_idmva, label='Signal')
    axs.hist(bins_idmva[:-1]+0.5*(bins_idmva[1:]-bins_idmva[:-1]),weights=bkg_tot_min_idmva, bins=bins_idmva, label='GJets')
    axs.hist(bins_idmva[:-1]+0.5*(bins_idmva[1:]-bins_idmva[:-1]),weights=bkg0_min_idmva, bins=bins_idmva, label='DiphotonJetsBox')
    axs.set_xlabel("Min IDMVA")
    axs.set_ylabel("Events")
    axs.legend(loc='best')
    fig.savefig('/afs/cern.ch/work/n/nschroed/ss_pyfit/CMSSW_10_2_14/src/diphoton-mva-optimizer/plots/dists_min_idmva.png')
    axs.set_ylim(10,120*max(bkg0_min_idmva.max(),bkg1_min_idmva.max()))
    plt.yscale("log")
    fig.savefig('/afs/cern.ch/work/n/nschroed/ss_pyfit/CMSSW_10_2_14/src/diphoton-mva-optimizer/plots/stacked_min_idmva_log.png')
    plt.close(fig)
    print("min_idmva")

    
    fig, axs = plt.subplots(nrows=1,ncols=1)
    axs.set_xlim(-1,1)
    axs.set_ylim(0.001,1.2*max(bkg0_min_idmva.max(),bkg1_min_idmva.max()))
    mids = bins_idmva[:-1]+0.5*(bins_idmva[1:]-bins_idmva[:-1])
    xloerr = mids - bins_idmva[:len(bins_idmva)-1]
    xhierr = bins_idmva[1:] - mids
    x_error = np.array([xloerr,xhierr])
    axs.errorbar(mids, bkg0_min_idmva, yerr=np.sqrt(bkg0_min_idmva), xerr=x_error, label="DiphotonJetsBox", color='dodgerblue', marker=",", mew=0, capsize=0, elinewidth=1,linestyle='')
    axs.errorbar(mids, bkg1_min_idmva, yerr=np.sqrt(bkg1_min_idmva), xerr=x_error, label="GJets", color='darkblue', marker=",", mew=0, capsize=0, elinewidth=1,linestyle='')
#    axs.errorbar(mids, sig_min_idmva*500, yerr=np.sqrt(sig_min_idmva*500), xerr=x_error, label="Signal*500", color='darkred', marker=",", mew=0, capsize=0, elinewidth=1,linestyle='')
    axs.set_xlabel("Min IDMVA")
    axs.set_ylabel("Events")
    axs.legend(loc='best')
    fig.savefig('/afs/cern.ch/work/n/nschroed/ss_pyfit/CMSSW_10_2_14/src/diphoton-mva-optimizer/plots/unstacked_min_idmva.png')
    axs.set_ylim(10,120*max(bkg0_min_idmva.max(),bkg1_min_idmva.max()))
    plt.yscale("log")
    fig.savefig('/afs/cern.ch/work/n/nschroed/ss_pyfit/CMSSW_10_2_14/src/diphoton-mva-optimizer/plots/unstacked_min_idmva_log.png')
    plt.close(fig)
    print('unstacked min_idmva')
    """


def plot(files):
    """
    plot a bunch of cross checks:
    R9 and pt distributions for signal and background
    """

    plot_stacked(files)
    plot_unstacked(files)

    return
