import pandas as pd
import numpy as np
import uproot as up


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

def extract_data(data_file):

    #columns needed for minimzation are hardcoded
    keep_cols = [
        'DiphotonMVA',
        'CMS_hgg_mass',
        'leadEta', 'subleadEta',
        'leadPt', 'subleadPt',
        'lead_R9', 'sublead_R9',
        'leadIDMVA', 'subleadIDMVA',
        'decorrSigmaM', 'weight']
    keep_cols_bkg = [
        'DiphotonMVA',
        'CMS_hgg_mass',
        'leadEta', 'subleadEta',
        'leadPt', 'subleadPt',
        'lead_R9', 'sublead_R9',
        'leadIDMVA', 'subleadIDMVA',
        'decorrSigmaM','weight',
        'gen_lead_pt', 'gen_sublead_pt']
    drop_cols = ['leadEta', 'subleadEta',
        'leadPt', 'subleadPt',
        'leadIDMVA', 'subleadIDMVA']

    #open and load the data
    config = open(args.inputFile,'r').readlines()
    config = [x.strip() for x in config]

    files = []

    for line in config:
        line_list = line.split('\t')
        print("[INFO] opening {} as dataframe".format(line_list[1]))
        df = pd.DataFrame()
        if line_list[0].find('Back') != -1 or line_list[0].find('back') != -1:
            df = up.open(line_list[1])[line_list[0]].pandas.df(keep_cols_bkg)
        else:
            df = up.open(line_list[1])[line_list[0]].pandas.df(keep_cols)
        #preselection
        files.append(df)

    data = pd.concat(df_collection_data)
    data.drop(drop_cols,axis=1,inplace=True)

    return df