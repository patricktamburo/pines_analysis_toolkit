import pathlib
import pdb 
from datetime import datetime
import numpy as np 
from astropy.io import fits

'''Authors: 
        Patrick Tamburo, Boston University, June 2020
   Purpose: 
        Idenfities the closest flat field in time to the image you want reduce, that's in the correct band.
   Inputs:
        flats_path (pathlib.Path): path to the local calibrations/flats directory. 
        header (astropy.io.fits.header.Header): the header of the image you want to reduce. 
    Outputs:
        master_flat (np.ndarray): 2D array containing the appropriate master flat data. 
        master_flat_name (str): the name of the chosen master flat file. 
    TODO:
        None
'''

def master_dark_stddev_chooser(dark_std_path, header):

    exptime = header['EXPTIME']
    obs_date = datetime.strptime(header['DATE-OBS'].split('T')[0].replace('-',''),'%Y%m%d')
    possible_dark_stds = [x for x in dark_std_path.glob('*'+str(exptime)+'*.fits')]
    if (len(possible_dark_stds)) == 0:
        print('ERROR: Could not find any suitable master_dark_stddev images to get read noise/dark current measurement for {}.'.format(header['FILENAME']))
        return
    else:
        possible_dark_std_dates = [datetime.strptime(i.name.split('_')[-1].split('.')[0],'%Y%m%d') for i in possible_dark_stds]
        dark_std_date_distances = [abs(possible_dark_std_dates[i]-obs_date) for i in range(len(possible_dark_std_dates))]
        dark_std_ind = np.where(np.array(dark_std_date_distances) == min(np.array(dark_std_date_distances)))[0][0]
        master_dark_std = fits.open(possible_dark_stds[dark_std_ind])[0].data
        master_dark_std_name = possible_dark_stds[dark_std_ind].name
        return master_dark_std