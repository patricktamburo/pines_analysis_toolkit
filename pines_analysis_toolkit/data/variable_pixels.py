import numpy as np
from pines_analysis_toolkit.utils import pines_dir_check, quick_plot as qp
import pdb
import natsort
from astropy.io import fits
from pathlib import Path
import pysftp
import os
from datetime import datetime
import time

'''Authors:
		Patrick Tamburo, Boston University, September 2020
	Purpose:
        Identifies variable pixels using a dark stddev image. 
	Inputs:
        date (str): the UT date during which the dome flat field data was obtained (i.e., '20200531')
        exptime (float): the exposure time of the dark images in question, in seconds
        upload (bool, optional): Whether or not to upload the master dark to pines.bu.edu. By default, set to False (so you will not attempt to upload!).
        sftp (pysftp.Connection, optional): Connection to the PINES server, needed if you want to upload variable pixel mask.
    Outputs:
		Saves vpm_exptime_s_date.fits to Calibrations/Variable Pixel Masks/
	TODO:
        None
'''
def variable_pixels(date, exptime, clip_lvl=5, upload=False, sftp=''):
    pines_path = pines_dir_check()
    dark_stddev_path = pines_path/('Calibrations/Darks/Master Darks Stddev/')
    all_dark_stddev_files = natsort.natsorted(list(Path(dark_stddev_path).rglob('*'+date+'.fits')))
    dark_stddev_file = ''
    for file in all_dark_stddev_files:
        if float(file.name.split('_')[3]) == exptime:
            dark_stddev_file = file

    if dark_stddev_file == '':
        raise RuntimeError('No dark stddev files found on disk with date '+date+' and exposure time '+str(exptime)+' seconds!')
    
    master_dark_stddev = fits.open(dark_stddev_file)[0].data
    shape = np.shape(master_dark_stddev)
    #Flag variable pixels, those that vary >= clip_lvl * the mean variation. Save a boolean mask of variable pixels to incorporate into the bad pixel mask. 

    variable_inds = np.where(master_dark_stddev >= clip_lvl*np.nanmean(master_dark_stddev))
    variable_mask = np.zeros((shape[0], shape[1]), dtype='int')
    variable_mask[variable_inds] = 1
    print('')
    print('Found {} variable pixels.'.format(len(variable_inds[0])))

    output_filename = 'vpm_'+str(exptime)+'_s_'+date+'.fits'
    output_path = pines_path/('Calibrations/Variable Pixel Masks/'+output_filename)

    #Add some header keywords detailing the master_dark creation process. 
    hdu = fits.PrimaryHDU(variable_mask)
    hdu.header['HIERARCH DATE CREATED'] = datetime.utcnow().strftime('%Y-%m-%d')+'T'+datetime.utcnow().strftime('%H:%M:%S')
    hdu.header['HIERARCH SIGMA CLIP LVL'] = clip_lvl
    hdu.header['HIERARCH CLIP VALUE'] = clip_lvl*np.nanmean(master_dark_stddev)

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
        sftp.chdir('data/calibrations/Variable Pixel Masks')
        upload_name = output_filename
        if upload_name in sftp.listdir():
            print('WARNING: This will overwrite {} in pines.bu.edu:data/calibrations/Variable Pixel Masks/'.format(upload_name))
            upload_check = input('Do you want to continue? y/n: ')
            if upload_check == 'y':
                sftp.put(output_path,upload_name)
                print('Uploaded to pines.bu.edu:data/calibrations/Variable Pixel Masks/!')
            else:
                print('Skipping upload!')
        else:
            sftp.put(output_path,upload_name)
            print('Uploaded {} to pines.bu.edu:data/calibrations/Variable Pixel Masks/!'.format(upload_name))