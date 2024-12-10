import pandas as pd
import numpy as np
import uproot3 as up
import logging

import python.plotters.plot as plot

def print_setup(cmd_line_args, input_file):
    """ prints the configuation of the program """
    cmd = ''
    for arg in cmd_line_args:
        if ' ' in arg:
            cmd += '"{}" '.format(arg)
        else:
            cmd +=  "{} ".format(arg)

    logging.info("#"*40)
    logging.info("[INFO] Welcome to the diphoton mva boundary optimizer")
    logging.info("[INFO] the command you ran was: {}".format(cmd))
    logging.info("[INFO] the specified config file is: {}".format(input_file))

def extract_data(args):
    """ loads config file into dataframes """
    #open args
    cfg_file = args.inputFile
    _kXcheckPlot = args.xcheck
    _kStackPlot = args.plot
    output_title = args.output
    lumi_scale = args.lumi
    lumi_scale_bkg = args.lumi_bkg
    bkg_scale = args.bkg_scale
    #open and load the data
    config = open(cfg_file, 'r').readlines()
    config = [x.strip() for x in config]

    keep_cols = ['CMS_hgg_mass', 'diphoton_transformedMva', 'diphoton_pt', 'weight']    
    keep_cols_sig = ['CMS_hgg_mass', 'diphoton_mva', 'diphoton_pt', 'weight']    

    bkg_files = []
    bkg_titles = []
    sig_files = []
    sig_titles = []

    for line in config:
        line_list = line.split(',')
        #line_list = line.split('\t')
        #print(len(line_list))
        #logging.info("[INFO] opening {} as dataframe".format(line_list[2]))
        #df = up.open(line_list[2])[line_list[1]].pandas.df(keep_cols)
        year_index = 1*('16' in line_list[3]) + 2*('17' in line_list[3]) + 3*('18' in line_list[3])
        if year_index == 0:
              raise ValueError("year index = {} does not correspond to 2016, 2017 or 2018, please indicate the year correctly in the legendEntry column in your config file".format(year_index))
        if line_list[0].find('bkg') != -1:
            logging.info("[INFO] opening {} as dataframe".format(line_list[2]))
            df = up.open(line_list[2])[line_list[1]].pandas.df(keep_cols)
            if ('QCD' in line_list[1]):
                df['weight'] = (1/40)* df['weight'].values*lumi_scale_bkg[year_index-1]
            else:
                df['weight'] = lumi_scale_bkg[year_index-1] * df['weight'].values
            bkg_files.append(df)
            bkg_titles.append(line_list[3])
        
        else:
            logging.info("[INFO] opening {} as dataframe".format(line_list[2]))
            df = up.open(line_list[2])[line_list[1]].pandas.df(keep_cols)
            df['weight'] = lumi_scale[year_index-1] * df['weight'].values
            sig_files.append(df)
            sig_titles.append(line_list[3])
        if _kXcheckPlot:
            plot.xcheck_plot(df, line_list[3])


    if _kStackPlot:
        sig = ([pd.concat(sig_files)], ["1000*signal"])
        bkg = (bkg_files, bkg_titles)
        plot.stack_plot(sig, bkg, output_title, lumi_scale, bounds=args.boundaries)

    df_bkg = pd.concat(bkg_files)
    df_sig = pd.concat(sig_files)
    #df_bkg.rename(columns={'diphoton_transformedMva': 'diphoton_mva'}, inplace=True)
    
    df_sig['is_signal'] = [True for i in range(len(df_sig))]
    df_bkg['is_signal'] = [False for i in range(len(df_bkg))]

    return df_sig, df_bkg
