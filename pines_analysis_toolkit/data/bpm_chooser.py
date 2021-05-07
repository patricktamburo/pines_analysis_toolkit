from datetime import datetime
import numpy as np
import pdb
from astropy.io import fits

'''Authors: 
        Patrick Tamburo, Boston University, September 2020
   Purpose: 
        Identifies the closest bpm in time to the image you want reduce, that has the correct exposure time and is in the right band.
   Inputs:
        bpm_path (pathlib.Path): path to the local Calibrations/Bad Pixel Masks/ directory. 
        header (astropy.io.fits.header.Header): the header of the image you want to reduce. 
    Outputs:
        master_bpm (np.ndarray): 2D array containing the appropriate master bad pixel mask. 
        master_bpm_name (str): the name of the chosen bad pixel mask file. 
    TODO:
        None
'''

def bpm_chooser(bpm_path, header):
    exptime = header['EXPTIME']
    band = header['FILTNME2']
    obs_date = datetime.strptime(header['DATE-OBS'].split('T')[0].replace('-',''),'%Y%m%d')
    possible_bpms = [x for x in (bpm_path.glob('*'+band+'_'+str(exptime)+'*.fits'))]

    if (len(possible_bpms) == 0): 
        print('Warning: Could not find any suitable bpms to reduce {}.'.format(header['FILENAME']))
        print('Expanding bpm search.')
        possible_bpms = [x for x in (bpm_path.glob('*'+band+'*.fits'))]
        if len(possible_bpms) == 0:
            print('ERROR: no suitable BPMs found inspect manually.')
            return
    possible_bpm_dates = [datetime.strptime(i.name.split('_')[-1].split('.')[0],'%Y%m%d') for i in possible_bpms]
    bpm_date_distances = [abs(possible_bpm_dates[i]-obs_date) for i in range(len(possible_bpm_dates))]
    bpm_ind = np.where(np.array(bpm_date_distances) == min(np.array(bpm_date_distances)))[0][0]
    master_bpm = fits.open(possible_bpms[bpm_ind])[0].data	
    master_bpm_name = possible_bpms[bpm_ind].name
    return master_bpm, master_bpm_name