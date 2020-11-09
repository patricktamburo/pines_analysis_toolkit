import pdb
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import glob
import pickle
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pathlib import Path
from datetime import datetime
import time
import os

'''Authors:
		Patrick Tamburo, Boston University, September 2020
	Purpose:
		Creates a combined bad pixel mask from Kokopelli, variable, hot, and dead pixel masks.
	Inputs:
		date (str): the date calibration data was taken 
        exptime (float): the exposure time of the darks, in seconds
        band (str): the band of the flats
		upload (bool, optional): whether or not to upload the bad pixel mask to pines.bu.edu. By default, False (so you won't try to upload!).
		sftp (pysftp.Connection, optional): the sftp connection to the pines server, required if you are going to upload bad pixel mask.
	Outputs:
		Writes bpm_band_exptime_s_date.fits to Calibrations/Bad Pixel Masks/
	TODO:
		None
'''

def bpm_maker(date, exptime, band, upload=False, sftp=''):
    pines_path = pines_dir_check()

    #Load in the different masks. 
    kokopelli_path = pines_path/('Calibrations/Kokopelli Mask/kokopelli_mask.fits')
    kokopelli_mask = (1-fits.open(kokopelli_path)[0].data).astype('int')[0:1024,:]

    variable_path = pines_path/('Calibrations/Variable Pixel Masks/vpm_'+str(exptime)+'_s_'+date+'.fits')
    variable_mask = fits.open(variable_path)[0].data

    hot_path = pines_path/('Calibrations/Hot Pixel Masks/hpm_'+str(exptime)+'_s_'+date+'.fits')
    hot_mask = fits.open(hot_path)[0].data

    dead_path = pines_path/('Calibrations/Dead Pixel Masks/dpm_'+band+'_'+date+'.fits')
    dead_mask = fits.open(dead_path)[0].data

    #Visualize all masks
    bpm = np.zeros(np.shape(dead_mask), dtype='int')
    bad_locs = np.where((kokopelli_mask == 1) | (variable_mask == 1) | (hot_mask == 1) | (dead_mask == 1))
    bpm[bad_locs] = 1
    
    num_bad = len(np.where(bpm ==1)[0])
    frac_bad = num_bad / 1024**2

    print('{} percent of the detector flagged as bad.'.format(np.round(frac_bad*100,1)))
    plt.ion()
    plt.imshow(bpm, origin='lower')


    output_filename = 'bpm_'+band+'_'+str(exptime)+'_s_'+date+'.fits'
    output_path = pines_path/('Calibrations/Bad Pixel Masks/'+output_filename)

    hdu = fits.PrimaryHDU(bpm)
    hdu.header['HIERARCH DATE CREATED'] = datetime.utcnow().strftime('%Y-%m-%d')+'T'+datetime.utcnow().strftime('%H:%M:%S')

    #Now save to a file on your local machine. 
    print('')
    print('Writing the file to '+output_filename)

    #Check to see if other files of this name exist.
    if os.path.exists(output_path):
        print('')
        print('WARNING: This will overwrite {}!'.format(output_path))
        dark_check = input('Do you want to continue? y/n: ')
        if dark_check == 'y':
            hdu.writeto(output_path,overwrite=True)
            print('Wrote to {}!'.format(output_path))
        else:
            print('Not overwriting!')
    else:
        hdu.writeto(output_path,overwrite=True)
    print('')

    #Upload the master dark to PINES server.
    if upload:
        print('Beginning upload process to pines.bu.edu...')
        print('Note, only PINES admins will be able to upload.')
        time.sleep(2)
        print('')
        sftp.chdir('/')
        sftp.chdir('data/calibrations/Bad Pixel Masks')
        upload_name = output_filename
        if upload_name in sftp.listdir():
            print('WARNING: This will overwrite {} in pines.bu.edu:data/calibrations/Bad Pixel Masks/'.format(upload_name))
            upload_check = input('Do you want to continue? y/n: ')
            if upload_check == 'y':
                sftp.put(output_path,upload_name)
                print('Uploaded to pines.bu.edu:data/calibrations/Bad Pixel Masks/!')
            else:
                print('Skipping upload!')
        else:
            sftp.put(output_path,upload_name)
            print('Uploaded {} to pines.bu.edu:data/calibrations/Bad Pixel Masks/!'.format(upload_name))
   