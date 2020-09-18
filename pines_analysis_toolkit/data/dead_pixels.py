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

def dead_pixels(date, band, clip_lvl=3, upload=False, sftp=''):
    pines_path = pines_dir_check()
    flat_path = pines_path/('Calibrations/Flats/Domeflats/'+band+'/Master Flats/')
    flat_file = natsort.natsorted(list(flat_path.rglob('*'+date+'.fits')))[0]
    master_flat = fits.open(flat_file)[0].data
    shape = np.shape(master_flat)

    print('')
    print('Flagging dead pixels')
    print('......')
    #Find hot pixels in the master dark. This uses the 8 pixels surrounding each pixel, ignoring the pixels on the edges of the detector.
    dead_pixel_mask = np.zeros((shape[0], shape[1]), dtype='int')
    num_flagged = 999 #Initialize
    total_flagged = 0
    iteration = 1
    max_iters = 5
    while (num_flagged != 0) or (iteration == max_iters): 
        num_flagged = 0
        pbar = ProgressBar()
        for xx in pbar(range(1, shape[0]-1)):
            for yy in range(1, shape[1]-1):
                if dead_pixel_mask[yy,xx] == False: #Only check pixels that aren't already flagged. 
                    pixel_val = master_flat[yy,xx]
                    neighbor_vals = np.delete(master_flat[yy-1:yy+2,xx-1:xx+2].ravel(),4) #This grabs the 8 surrounding pixels, ignoring the target pixel. 
                    if pixel_val < np.nanmean(neighbor_vals) - clip_lvl*np.nanstd(neighbor_vals): #Flag pixel as hot if it is more than clip_lvl sigma higher than the mean of its neighbors. 
                        dead_pixel_mask[yy,xx] = 1
                        num_flagged += 1
        print('Iteration {}: {} new dead pixels identified.'.format(iteration, num_flagged))
        iteration += 1
        total_flagged += num_flagged

    print('')
    print('Found {} dead pixels.'.format(total_flagged))

    output_filename = 'dpm_'+band+'_'+date+'.fits'
    output_path = pines_path/('Calibrations/Dead Pixel Masks/'+output_filename)

    hdu = fits.PrimaryHDU(dead_pixel_mask)
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
        sftp.chdir('data/calibrations/Dead Pixel Masks')
        upload_name = output_filename
        if upload_name in sftp.listdir():
            print('WARNING: This will overwrite {} in pines.bu.edu:data/calibrations/Dead Pixel Masks/'.format(upload_name))
            upload_check = input('Do you want to continue? y/n: ')
            if upload_check == 'y':
                sftp.put(output_path,upload_name)
                print('Uploaded to pines.bu.edu:data/calibrations/Dead Pixel Masks/!')
            else:
                print('Skipping upload!')
        else:
            sftp.put(output_path,upload_name)
            print('Uploaded {} to pines.bu.edu:data/calibrations/Dead Pixel Masks/!'.format(upload_name))