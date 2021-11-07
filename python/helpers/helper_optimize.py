import pandas as pd
import numpy as np
import uproot3 as up

import python.plotters.plot as plot

def print_setup(cmd_line_args, input_file):
    cmd = ''
    for arg in cmd_line_args:
        if ' ' in arg:
            cmd += '"{}" '.format(arg)
        else:
            cmd +=  "{} ".format(arg)

    print("#"*40)
    print("[INFO] Welcome to the diphoton mva boundary optimizer")
    print("[INFO] the command you ran was: {}".format(cmd))
    print("[INFO] the specified config file is: {}".format(input_file))

def extract_data(cfg_file, _kXcheckPlot, _kStackPlot, output_title):

    #open and load the data
    config = open(cfg_file, 'r').readlines()
    config = [x.strip() for x in config]

    keep_cols = ['CMS_hgg_mass', 'diphoton_mva', 'weight']

    bkg_files = []
    bkg_titles = []
    sig_files = []
    sig_titles = []

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
            sig_files.append(df)
            sig_titles.append(line_list[3])
        if _kXcheckPlot:
            plot.xcheck_plot(df, line_list[3])

    if _kStackPlot:
        sig = (sig_files, sig_titles)
        bkg = (bkg_files, bkg_titles)
        plot.stack_plot(sig, bkg, output_title)

    df_bkg = pd.concat(bkg_files)
    df_sig = pd.concat(sig_files)

    cols_bkg = list(df_bkg.columns)
    cols_sig = list(df_sig.columns)

    drop_cols = [col for col in cols_bkg + cols_sig if (col not in cols_bkg) or (col not in cols_sig)]
    drop_cols_sig = [col for col in drop_cols if col in cols_sig]
    drop_cols_bkg = [col for col in drop_cols if col in cols_bkg]

    if len(drop_cols_bkg) > 0:
        df_bkg.drop(drop_cols_bkg, axis=1, inplace=True)
    if len(drop_cols_sig) > 0:
        df_sig.drop(drop_cols_sig, axis=1, inplace=True)

    return df_sig, df_bkg
