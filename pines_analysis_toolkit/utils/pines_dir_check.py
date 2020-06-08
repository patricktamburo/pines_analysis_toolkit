from pathlib import Path
import os
import pdb

def pines_dir_check():
    home_dir = Path(os.path.expanduser('~/Documents/'))
    default_pines_dir = os.path.join(home_dir,'PINES_analysis_toolkit/')
    if os.path.exists(default_pines_dir):
        return default_pines_dir
    else:
        print('ERROR...I have not set this up to work for directories other than ~/Documents/PINES_analysis_toolkit.')
        #TODO: Config file for non-Documents path? 
        return