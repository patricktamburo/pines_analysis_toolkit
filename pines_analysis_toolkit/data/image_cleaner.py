import pines_analysis_toolkit as pat
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from astropy.stats import sigma_clipped_stats
import pdb 
import copy
import time
import pickle
from scipy.stats import sigmaclip
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.data.bg_2d import bg_2d
import natsort
from photutils import CircularAperture, CircularAnnulus
import os 
'''Authors:
		Patrick Tamburo, Boston University, August 2020
	Purpose:
		Identifies bad pixels in PINES images by examing each pixel's value in relation to the mean and standard deviation of its neighbors. 
	Inputs:
		target_name (str): the target's full 2MASS name, i.e. '2MASS J01234567+0123456' 
        centroided_sources (pandas dataframe): source centroid positions in every image you want to measure, output from pat.photometry.centroider. 
        r_out (int, optional): The outer annulus radius that will be used in photometry. pixels falling within this radius of each source will be checked to see if they're bad. 
        box_l (int, optional): the side length of the box in pixels, with the target pixel at the center. NOTE: Must be an odd integer!
        sigma (float, optional): the number of standard deviations that a pixel can be discrepant from the mean of its neighbors before being flagged as bad. 
	Outputs:
		Saves bad-pixel-flagged images to Object's bp_reduced/ directory. Bad pixels are replaced with NaNs. 
    TODO: 
        Only need to worry about flagging pixels in the vicinity of our sources, will speed things up a lot. 
'''
def image_cleaner(target, centroided_sources, r_out=30, box_l=5, sigma=4.0):
    if box_l % 2 != 1:
        raise ValueError('box_w must be an odd integer!')

    #Set up paths. 
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    red_path = pines_path/('Objects/'+short_name+'/reduced/')
    red_files = np.array(natsort.natsorted([x for x in red_path.glob('*.fits')]))
    output_path = pines_path/('Objects/'+short_name+'/bp_bg_reduced/')

    #Get source names. 
    centroided_sources.columns = centroided_sources.columns.str.strip()
    source_names = natsort.natsorted(list(set([i[0:-2].replace('X','').replace('Y','').rstrip().lstrip() for i in centroided_sources.keys()])))
    if len(source_names) != len(centroided_sources.keys())/2:
        print('ERROR: something going wrong grabbing source names. ')
        return
    if len(centroided_sources) != len(red_files):
        print('ERROR: number of centroided images does not match number of reduced images.')
        return 
    
    #Load in the bad pixel mask created from biases/flats. 
    bpm_path = '/Users/tamburo/Documents/PINES_analysis_toolkit/Calibrations/Bad Pixel Masks/bpm.p'
    bpm = pickle.load(open(bpm_path, 'rb'))

    #Loop over all reduced files and identify bad pixels near the sources in each. 
    for i in range(len(red_files)):
        file = red_files[i]
        #Read in the image/header.
        data = fits.open(file)[0].data
        header = fits.open(file)[0].header
        target_filename = output_path/(file.name.split('.fits')[0]+'_bp_bg.fits')
        if not os.path.exists(target_filename):
            #Use the bad pixel mask to preemptively NaN pixels in the image. 
            data[np.where(bpm == 0)] = np.nan

            print('{}, image {} of {}.'.format(file.name, i+1, len(red_files)))

            #Loop over all sources.  
            for source in source_names:
                source_x = int(np.round(centroided_sources[source+' X'][i]))
                source_y = int(np.round(centroided_sources[source+' Y'][i]))

                y_len, x_len = r_out*2, r_out*2 #The side lenght of the box around the source to check for bad pixels. 
                x_start, y_start = int(source_x - x_len/2), int(source_y - y_len/2)
                x_end, y_end = int(source_x + x_len/2), int(source_y + y_len/2)

                #Handle reference stars near the edges. 
                if x_end > 1023:
                    x_end = 1023
                if y_end > 1023:
                    y_end = 1023
                if x_start < 0:
                    x_start = 0
                if y_start < 0:
                    y_start = 0

                #Loop pixels surrounding the source position, cut out boxes around them, and check if they are more than sigma standard deviations away from their neighbors. 
                num_bad = 9999 #Initialize
                loop = 1

                #Loop over all pixels in the image, correcting them if it's >sigma away from the mean of its neighbors. 
                #Iterate until no new bad pixels are flagged.
                while num_bad != 0:
                    num_bad = 0 
                    for x in range(x_start, x_end+1):
                        for y in range(y_start, y_end+1): 
                            if np.sqrt((x-source_x)**2 + (y-source_y)**2) < r_out: #Only check pixels within r_out of the target_position. 
                                cutout = []
                                target_pixel_val = data[y, x]

                                #Get the start/stop pixels for a box of side length box_l surrounding the target pixel.
                                pix_box_y_start = max(0, y-int(box_l/2))
                                pix_box_y_end = min(1024, y+int(box_l/2)+1)
                                pix_box_x_start = max(0, x-int(box_l/2))
                                pix_box_x_end = min(1024, x+int(box_l/2)+1)
                                    
                                for y_cut in range(pix_box_y_start, pix_box_y_end):
                                    for x_cut in range(pix_box_x_start, pix_box_x_end):
                                        if not ((y_cut == y) and (x_cut == x)) and not np.isnan(data[y_cut, x_cut]):
                                            #Append the neighboring pixel values if they are not NaNs. 
                                            cutout.append(data[y_cut,x_cut])
                                
                                #Do a sigmaclip on the cutout, in case it includes any bad pixels not yet flagged!
                                vals, lo, hi = sigmaclip(cutout, low=3., high=3.)

                                if abs(target_pixel_val - np.mean(vals)) / np.std(vals) > sigma:
                                    #If the target pixel value deviates by more than sigma above the mean of the neighboring pixels, make it a NaN. 
                                    num_bad += 1
                                    #If the bad pixel is >5 pixels away from the source centroid, replace it with the median of its neighbors. 
                                    if np.sqrt((x-source_x)**2 + (y-source_y)**2) > 5: 
                                        data[y,x] = np.median(vals)
                                    else:
                                        data[y,x] = np.nan
                    loop += 1

            #Remove large background trends. 
            #data = bg_2d(data)

            #Save the image. 
            target_filename = output_path/(file.name.split('.fits')[0]+'_bp_bg.fits')
            fits.writeto(target_filename, data, header, overwrite=True)
        else:
            print(file.name, ' already exists, skipping. ')
