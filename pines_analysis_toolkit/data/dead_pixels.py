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
import matplotlib.pyplot as plt
from scipy.stats import sigmaclip
from astropy.stats import histogram

'''Authors:
		Patrick Tamburo, Boston University, September 2020
	Purpose:
		Creates a dead pixel mask from master flat. 
	Inputs:
		date (str): the date calibration data was taken 
        band (str): the band of the flats
        clip_lvl (float, optional): the sigma level below which to flag dead pixels. 
		upload (bool, optional): whether or not to upload the dead pixel mask to pines.bu.edu. By default, False (so you won't try to upload!).
		sftp (pysftp.Connection, optional): the sftp connection to the pines server, required if you are going to upload dead pixel mask.
	Outputs:
		Writes dpm_band_date.fits to Calibrations/Dead Pixel Masks/
	TODO:
		None
'''

def dead_pixels(date, band, clip_lvl=4.5, box_l=3, upload=False, sftp=''):
    pines_path = pines_dir_check()
    flat_path = pines_path/('Calibrations/Flats/Domeflats/'+band+'/Master Flats/')
    flat_file = natsort.natsorted(list(flat_path.rglob('*'+date+'.fits')))[0]
    master_flat = fits.open(flat_file)[0].data
    shape = np.shape(master_flat)

    # #Can plot a histogram of master_flat pixel values. 
    # his = histogram(master_flat, bins=250)
    # plt.figure(figsize=(12,5))
    # plt.bar(his[1][0:-1], his[0], width=0.005) 
    # plt.yscale('log')
    # plt.xlabel('Pixel Value', fontsize=14)
    # plt.ylabel('N$_{pix}$', fontsize=14)
    # plt.title('Distribution of Pixel Values in Master Flat '+band+' '+date, fontsize=16)
    # pdb.set_trace()

    #Incorporate bad pixel information from the Kokopelli pixel mask. 
    kokopelli_path = pines_path/('Calibrations/Kokopelli Mask/kokopelli_mask.fits')
    kokopelli_mask = (1-fits.open(kokopelli_path)[0].data).astype('int')[0:1024,:]
    master_flat[np.where(kokopelli_mask == 1)] = np.nan
    
    print('')
    print('Flagging dead pixels.')
    
    #Find hot pixels in the master dark. This uses the 8 pixels surrounding each pixel, ignoring the pixels on the edges of the detector.
    #This will iterate until no new dead pixels are found. 
    dead_pixel_mask = np.zeros((shape[0], shape[1]), dtype='int')
    total_flagged = 0

    box_ls = [13, 11, 9, 7, 5]
    clip_lvls = [4.5, 4.5, 4.5, 4, 4]
    iteration = 1
    for i in range(len(box_ls)):
        box_l = box_ls[i]
        clip_lvl = clip_lvls[i]
        print('Box size = {} x {}.'.format(box_l, box_l))
        print('Sigma clipping level = {}.'.format(clip_lvl))
        print('......')
        num_flagged = 999 #Initialize
        box_minus = int(box_l/2)
        box_plus = int(box_l/2) + 1
        targ_pix_ind = int(box_l**2/2)

        while (num_flagged > 1): 
            num_flagged = 0
            pbar = ProgressBar()
            for xx in pbar(range(int(box_l/2), shape[0]-int(box_l/2))):
                for yy in range(int(box_l/2), shape[1]-int(box_l/2)):
                    if (not np.isnan(master_flat[yy,xx])): #Only check pixels that aren't already flagged. 
                        box_2d = master_flat[yy-box_minus:yy+box_plus,xx-box_minus:xx+box_plus] #2d cutout surrounding the target pixel.
                        box_1d = box_2d[~np.isnan(box_2d)].ravel() #Unraveled box, ignoring any NaNs. 
                        neighbor_vals = np.delete(box_1d, targ_pix_ind) #Remove the target pixel from the 1d array. 

                        #neighbor_vals = sigmaclip(n, low=3.5, high=3.5)[0] #Clip away other bad pixels in the box to avoid biasing the mean/standard deviation.                   
                        if master_flat[yy,xx] < (np.mean(neighbor_vals) - clip_lvl*np.std(neighbor_vals)): #Flag pixel as hot if it is more than clip_lvl sigma higher than the mean of its neighbors.
                            dead_pixel_mask[yy,xx] = 1
                            master_flat[yy,xx] = np.nan #Set this pixel to a NaN so that it's ignored on subsequent iterations. 
                            num_flagged += 1
            iteration += 1
            total_flagged += num_flagged
            print('Iteration {}: {} new dead pixels identified, {} dead pixels total.'.format(iteration, num_flagged, total_flagged))
            print('')

    print('')
    print('Found {} dead pixels.'.format(total_flagged))

    output_filename = 'dpm_'+band+'_'+date+'.fits'
    output_path = pines_path/('Calibrations/Dead Pixel Masks/'+output_filename)

    hdu = fits.PrimaryHDU(dead_pixel_mask)
    hdu.header['HIERARCH DATE CREATED'] = datetime.utcnow().strftime('%Y-%m-%d')+'T'+datetime.utcnow().strftime('%H:%M:%S')
    hdu.header['HIERARCH SIGMA CLIP LVL'] = clip_lvl

    #Now save to a file on your local machine. 
    print('')
    print('Writing the file to '+output_filename)

    #Check to see if other files of this name exist.
    # if os.path.exists(output_path):
    #     print('')
    #     print('WARNING: This will overwrite {}!'.format(output_path))
    #     dark_check = input('Do you want to continue? y/n: ')
    #     if dark_check == 'y':
    #         hdu.writeto(output_path,overwrite=True)
    #         print('Wrote to {}!'.format(output_path))
    #     else:
    #         print('Not overwriting!')
    # else:
    hdu.writeto(output_path,overwrite=True)
    print('Wrote to {}!'.format(output_path))
    print('')

    #Upload the master dark to PINES server.
    if upload:
        print('Beginning upload process to pines.bu.edu...')
        print('Note, only PINES admins will be able to upload.')
        print('')
        sftp.chdir('/')
        sftp.chdir('data/calibrations/Dead Pixel Masks')
        upload_name = output_filename
        # if upload_name in sftp.listdir():
        #     print('WARNING: This will overwrite {} in pines.bu.edu:data/calibrations/Dead Pixel Masks/'.format(upload_name))
        #     upload_check = input('Do you want to continue? y/n: ')
        #     if upload_check == 'y':
        #         sftp.put(output_path,upload_name)
        #         print('Uploaded to pines.bu.edu:data/calibrations/Dead Pixel Masks/!')
        #     else:
        #         print('Skipping upload!')
        # else:
        sftp.put(output_path,upload_name)
        print('Uploaded {} to pines.bu.edu:data/calibrations/Dead Pixel Masks/!'.format(upload_name))