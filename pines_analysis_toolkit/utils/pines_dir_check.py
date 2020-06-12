from pathlib import Path
import os
import pdb

'''Authors:
		Patrick Tamburo, Boston University, June 2020
	Purpose:
        Finds the local PINES_analysis_toolkit directory
	Inputs:
		None
    Outputs:
		default_pines_dir (pathlib.PosixPath): The path to the PINES_analysis_toolkit directory
	TODO:
        Get version working in case that the directory is NOT in ~/Documents. May have to use a config.txt file. 
'''

def pines_dir_check():
    home_dir = Path(os.path.expanduser('~/Documents/'))
    default_pines_dir = Path(os.path.join(home_dir,'PINES_analysis_toolkit/'))
    if os.path.exists(default_pines_dir):
        return default_pines_dir
    else:
        print('ERROR...I have not set this up to work for directories other than ~/Documents/PINES_analysis_toolkit.')
        #TODO: Config file for non-Documents path? 
        return