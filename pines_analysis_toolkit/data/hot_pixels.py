import numpy as np
from progressbar import ProgressBar
from pines_analysis_toolkit.utils import pines_dir_check
import natsort
from pathlib import Path
import pdb
from astropy.io import fits
from datetime import datetime 
import time 
import os 
from pines_analysis_toolkit.utils.quick_plot import quick_plot as qp

'''Authors:
		Patrick Tamburo, Boston University, September 2020
	Purpose:
		Creates a hot pixel mask from master dark. 
	Inputs:
		date (str): the date calibration data was taken 
        exptime (float): the exposure time of the darks in seconds.
        clip_lvl (float, optional): the sigma level above which to flag hot pixels. 
		upload (bool, optional): whether or not to upload the hot pixel mask to pines.bu.edu. By default, False (so you won't try to upload!).
        box_l (int, optional): The side length of the box of neighboring pixels to consider when flagging hot pixels. 
		sftp (pysftp.Connection, optional): the sftp connection to the pines server, required if you are going to upload hot pixel mask.
	Outputs:
		Writes hpm_exptime_s_date.fits to Calibrations/Dead Pixel Masks/
	TODO:
		None
'''

def hot_pixels(date, exptime, clip_lvl=6, box_l=3, saturation=4000, upload=False, sftp=''):

    if box_l % 2 == 0:
        raise RuntimeError('box_l must be odd!')

    pines_path = pines_dir_check()
    darks_path = pines_path/('Calibrations/Darks/Master Darks/')
    all_dark_stddev_files = natsort.natsorted(list(Path(darks_path).rglob('*'+date+'.fits')))
    dark_file = ''
    for file in all_dark_stddev_files:
        if float(file.name.split('_')[2]) == exptime:
            dark_file = file
    
    if dark_file == '':
        raise RuntimeError('No dark stddev files found on disk with date '+date+' and exposure time '+str(exptime)+' seconds!')
    
    master_dark = fits.open(dark_file)[0].data
    shape = np.shape(master_dark)
    hot_pixel_mask = np.zeros((shape[0], shape[1]), dtype='int')
    #Mask/nan any pixel that is > saturation in the master_dark
    hot_pixel_mask[np.where(master_dark > saturation)] = 1
    master_dark[np.where(master_dark > saturation)] = np.nan
    
    #Incorporate information from the variable/Kokopelli pixel masks. 
    variable_path = pines_path/('Calibrations/Variable Pixel Masks/vpm_'+str(exptime)+'_s_'+date+'.fits')
    variable_mask = fits.open(variable_path)[0].data

    kokopelli_path = pines_path/('Calibrations/Kokopelli Mask/kokopelli_mask.fits')
    kokopelli_mask = (1-fits.open(kokopelli_path)[0].data).astype('int')[0:1024,:]

    master_dark[np.where(kokopelli_mask == 1)] = np.nan
    master_dark[np.where(variable_mask == 1)] = np.nan

    print('')
    print('Flagging hot pixels')
    print('......')
    #Find hot pixels in the master dark. This uses the 8 pixels surrounding each pixel, ignoring the pixels on the edges of the detector.
    #This will iterate until it finds no more hot pixels. 
    num_flagged = 999 #Initialize
    total_flagged = 0
    iteration = 1
    while (num_flagged != 0): 
        num_flagged = 0
        pbar = ProgressBar()
        for xx in pbar(range(int(box_l/2), shape[0]-int(box_l/2))):
            for yy in range(int(box_l/2), shape[1]-int(box_l/2)):
                if (hot_pixel_mask[yy,xx] == 0) and not np.isnan(master_dark[yy,xx]):
                    pixel_val = master_dark[yy,xx]
                    neighbor_vals = np.delete(master_dark[yy-int(box_l/2):yy+int(box_l/2)+1,xx-int(box_l/2):xx+int(box_l/2)+1].ravel(),4) #This grabs the surrounding pixels, ignoring the target pixel. 
                    if pixel_val > np.nanmean(neighbor_vals) + clip_lvl*np.nanstd(neighbor_vals): #Flag pixel as hot if it is more than clip_lvl sigma higher than the mean of its neighbors. 
                        hot_pixel_mask[yy,xx] = 1
                        master_dark[yy,xx] = np.nan #Set this pixel to a NaN so that it's ignored on subsequent iterations. 
                        num_flagged += 1
        print('Iteration {}: {} new hot pixels identified.'.format(iteration, num_flagged))
        iteration += 1
        total_flagged += num_flagged

    print('')
    print('Found {} hot pixels.'.format(total_flagged))

    output_filename = 'hpm_'+str(exptime)+'_s_'+date+'.fits'
    output_path = pines_path/('Calibrations/Hot Pixel Masks/'+output_filename)

    #Add some header keywords detailing the master_dark creation process. 
    hdu = fits.PrimaryHDU(hot_pixel_mask)
    hdu.header['HIERARCH DATE CREATED'] = datetime.utcnow().strftime('%Y-%m-%d')+'T'+datetime.utcnow().strftime('%H:%M:%S')
    hdu.header['HIERARCH SIGMA CLIP LVL'] = clip_lvl

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
        sftp.chdir('data/calibrations/Hot Pixel Masks')
        upload_name = output_filename
        if upload_name in sftp.listdir():
            print('WARNING: This will overwrite {} in pines.bu.edu:data/calibrations/Hot Pixel Masks/'.format(upload_name))
            upload_check = input('Do you want to continue? y/n: ')
            if upload_check == 'y':
                sftp.put(output_path,upload_name)
                print('Uploaded to pines.bu.edu:data/calibrations/Hot Pixel Masks/!')
            else:
                print('Skipping upload!')
        else:
            sftp.put(output_path,upload_name)
            print('Uploaded {} to pines.bu.edu:data/calibrations/Hot Pixel Masks/!'.format(upload_name))