import pdb
def get_source_names(df):
    ''' 
        Author: 
            Patrick Tamburo, Boston University, February 2021.
        Purpose:
            Grabs sources names from a PINES centroided sources dataframe.
        Inputs: 
            df (pandas DataFrame): dataframe output from centroider()
        Outputs:
            source_names (list): List of source names. 
    '''

    df_keys = df.keys().str.strip().str.replace(' Corrected Flux Error', '').str.replace(' Image X','').str.replace(' Image Y','').str.replace(' Cutout X','').str.replace(' Cutout Y','').str.replace(' Centroid Warning', '').str.replace(' Corrected Flux','').str.replace(' ALC Weight', '')
    source_names = []
    for i in range(len(df_keys)):
        df_key = df_keys[i]
        if (df_key != 'Filename') and (df_key != 'Seeing') and (df_key != 'Airmass') and (df_key != 'Time BJD TDB') and (df_key != 'Sigma Clip Flag') and (df_key not in source_names):
            source_names.append(df_key)
    return source_names