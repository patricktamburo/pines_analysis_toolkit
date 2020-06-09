import pandas

'''Authors:
		Patrick Tamburo, Boston University, June 2020
	Purpose:
        Reads in a pines YYYYMMDD_log.txt file and outputs a pandas dataframe of the log. 
	Inputs:
		path (str): the path to the log on your local machine
    Outputs:
		None
	TODO:
		None
'''

def pines_log_reader(path):
    df = pandas.read_csv(path)

    #Remove trailing/leading spaces from column names
    df.columns = df.columns.str.lstrip()
    df.columns = df.columns.str.rstrip()

    #Remove our header comment idicator in the first column if it's there.
    if '#' in df.columns[0]:
        df.rename(columns={df.columns[0]:df.columns[0].replace('#','')}, inplace=True)
    
    #Remove trailing and leading spaces from log entries. 
    for key in df.keys():
        try:
            df[key] = df[key].str.strip()
        except:
            continue
        
    return df