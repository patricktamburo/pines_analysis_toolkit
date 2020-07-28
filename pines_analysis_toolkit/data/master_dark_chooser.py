from datetime import datetime
import numpy as np
from astropy.io import fits
import pathlib
import pdb

'''Authors: 
        Patrick Tamburo, Boston University, June 2020
   Purpose: 
        Idenfities the closest dark in time to the image you want reduce, that has the correct exposure time.
   Inputs:
        dark_path (pathlib.Path): path to the local calibrations/darks directory. 
        header (astropy.io.fits.header.Header): the header of the image you want to reduce. 
    Outputs:
        master_dark (np.ndarray): 2D array containing the appropriate master dark data. 
        master_dark_name (str): the name of the chosen master dark file. 
    TODO:
        None
'''

def master_dark_chooser(dark_path, header):
    exptime = header['EXPTIME']
    obs_date = datetime.strptime(header['DATE-OBS'].split('T')[0].replace('-',''),'%Y%m%d')
    possible_darks = [x for x in (dark_path/'Master Darks').glob('*'+str(exptime)+'*.fits')]
    if (len(possible_darks) == 0): 
        print('ERROR: Could not find any suitable darks to reduce {}.'.format(header['FILENAME']))
        return
    else: 
        possible_dark_dates = [datetime.strptime(i.name.split('_')[-1].split('.')[0],'%Y%m%d') for i in possible_darks]
        dark_date_distances = [abs(possible_dark_dates[i]-obs_date) for i in range(len(possible_dark_dates))]
        dark_ind = np.where(np.array(dark_date_distances) == min(np.array(dark_date_distances)))[0][0]
        master_dark = fits.open(possible_darks[dark_ind])[0].data	
        master_dark_name = possible_darks[dark_ind].name
        return master_dark, master_dark_name