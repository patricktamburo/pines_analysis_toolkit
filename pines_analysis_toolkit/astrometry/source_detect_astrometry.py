import pines_analysis_toolkit as pat 
from glob import glob 
from natsort import natsorted 
import numpy as np
from astropy.convolution import interpolate_replace_nans, Gaussian2DKernel
from astropy.io import fits 
import os 
import time
from pathlib import Path

def source_detect_astrometry(target, api_key, filename, download_data=False):
    '''
        Authors:
            Pat Tamburo, Boston University, May 2021
        Purpose: 
            Uploads a single reduced image for a target to astrometry.net, downloads solution image, and updates the image header with the astrometry.net wcs. 
        Inputs:
            target (str): The target's long 2MASS name (e.g., '2MASS J09161504+2139512')
            api_key (str): The api key of your account on astrometry.net. 
            filename (pathlib.Path): Path to the source_detect_image.
            download_data (bool, optional): Whether or not to first download reduced data for the target from the PINES server. 
    '''

    kernel = Gaussian2DKernel(x_stddev=0.25)

    #If the header does not have a HISTORY keyword (which is added by astrometry.net), process it. 
    header = fits.open(filename)[0].header

    if 'HISTORY' not in header:
        print('Uploading {} to astrometry.net.'.format(filename.name))
        #Read in the image data, interpolate NaNs, and save to a temporary fits file. 
        #Astrometry.net does not work with NaNs in images. 
        original_image = fits.open(filename)[0].data
        image = interpolate_replace_nans(original_image, kernel=kernel)
        temp_filename = filename.parent/(filename.name.split('.fits')[0]+'_temp.fits')
        hdu = fits.PrimaryHDU(image, header=header)
        hdu.writeto(temp_filename, overwrite=True)

        #Upload the image with interpolated NaNs to astrometry.net. It will solve and download to a new image. 
        pat.astrometry.upload_to_astrometry.upload_file(api_key, temp_filename, header)

        #Try to donwload the solved image and open it. 
        try:
            #Grab the header of the astrometry.net solution image, and the original image data. 
            astrometry_image_path = filename.parent/(temp_filename.name.split('.fits')[0]+'_new_image.fits')
            wcs_header = fits.open(astrometry_image_path)[0].header
            wcs_hdu = fits.PrimaryHDU(original_image, header=wcs_header)

            #Save the original image data with the new wcs header. 
            output_filename = filename
            wcs_hdu.writeto(output_filename, overwrite=True)

        #If the try clause didn't work, that's because the processing on astrometry.net failed. 
        except:
            raise RuntimeError('Astrometry solution failed! Use a different source detect image for choosing reference stars.')            

        #Delete temporary files. 
        os.remove(temp_filename)
        os.remove(astrometry_image_path)
        time.sleep(1)
        print('')

    #If the header DOES have a HISTORY keyword, skip it, it has already been processed. 
    else:
        print('Astrometric processing already complete for {}, skipping.'.format(filename.name))
        print('')