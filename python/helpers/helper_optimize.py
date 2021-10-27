import pandas as pd
import numpy as np
import uproot as up

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

def extract_data(cfg_file, _kPlot):

    #open and load the data
    config = open(cfg_file,'r').readlines()
    config = [x.strip() for x in config]

    bkg_files = []
    sig_files = []

    for line in config:
        line_list = line.split('\t')
        print("[INFO] opening {} as dataframe".format(line_list[1]))
        df = pd.DataFrame()
        if line_list[0].find('bkg') != -1:
            df = up.open(line_list[2])[line_list[1]].pandas.df()
            bkg_files.append(df)
        else:
            df = up.open(line_list[2])[line_list[1]].pandas.df()
            sig_files.append(df)
        if _kPlot:
            plot.plot(df, line_list[3])
        #preselection
    
    df_bkg = pd.concat(bkg_files)
    df_sig = pd.concat(sig_files)

    cols_bkg = list(df_bkg.columns)
    cols_sig = list(df_sig.columns)
    print(cols_bkg)
    print(cols_sig)
    drop_cols = [col for col in cols_bkg+cols_sig if (col not in cols_bkg) or (col not in cols_sig)]
    drop_cols_sig = [col for col in drop_cols if col in cols_sig]
    drop_cols_bkg = [col for col in drop_cols if col in cols_bkg]

    print(drop_cols_sig)
    print(drop_cols_bkg)

    if len(drop_cols_bkg) > 0:
        df_bkg.drop(drop_cols_bkg, axis=1, inplace=True)
    if len(drop_cols_sig) > 0:
        df_sig.drop(drop_cols_sig, axis=1, inplace=True)

    

    return df_sig, df_bkg
