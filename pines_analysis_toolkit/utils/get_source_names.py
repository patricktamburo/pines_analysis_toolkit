def get_source_names(centroided_sources):
    ''' 
        Author: 
            Patrick Tamburo, Boston University, February 2021.
        Purpose:
            Grabs sources names from a PINES centroided sources dataframe.
        Inputs: 
            centroided_souces (pandas DataFrame): dataframe output from centroider()
        Outputs:
            source_names (list): List of source names. 
    '''
    df_keys = centroided_sources.keys().str.strip().str.replace(' Image X','').str.replace(' Image Y','').str.replace(' Cutout X','').str.replace(' Cutout Y','').str.replace(' Centroid Warning', '')
    source_names = []
    for i in range(len(df_keys)):
        df_key = df_keys[i]
        if (df_key != 'Filename') and (df_key != 'Seeing') and (df_key != 'Airmass') and (df_key != 'Time (JD UTC)') and (df_key not in source_names):
            source_names.append(df_key)
    return source_names