import pines_analysis_toolkit as pat 
from glob import glob 
from natsort import natsorted 
import numpy as np
import pdb 
from astropy.convolution import interpolate_replace_nans, Gaussian2DKernel
from astropy.io import fits 
import os 

target = '2MASS J09130320+1841501'
api_key = 'vqrhxpyypnpmqbzo'


pines_path = pat.utils.pines_dir_check()
short_name = pat.utils.short_name_creator(target)
reduced_path = pines_path/('Objects/'+short_name+'/reduced/')
reduced_files = np.array(natsorted(glob(str(reduced_path)+'/*.fits')))

kernel = Gaussian2DKernel(x_stddev=0.25)


for i in range(len(reduced_files)):
    filename = reduced_files[i]

    #Read in the file, interpolate NaNs, and save to a temporary fits file. 
    #Astrometry.net does not work with NaNs in images. 
    image = interpolate_replace_nans(fits.open(filename)[0].data, kernel=kernel)
    header = fits.open(filename)[0].header
    temp_filename = filename.split('.fits')[0]+'_temp.fits'
    hdu = fits.PrimaryHDU(image, header=header)
    hdu.writeto(temp_filename, overwrite=True)

    #Upload the image with interpolated NaNs to astrometry.net. It will solve and download to a new image. 
    pat.astrometry.upload_to_astrometry.upload_file(api_key, temp_filename)

    #Grab the header of the astrometry.net solution image, and the original image data. 
    astrometry_image_path = temp_filename.split('.fits')[0]+'_new_image.fits'
    wcs_header = fits.open(astrometry_image_path)[0].header
    original_image = fits.open(filename)[0].data
    wcs_hdu = fits.PrimaryHDU(original_image, header=wcs_header)

    #Save the original image data with the new wcs header. 
    output_filename = filename.split('.fits')[0]+'_wcs.fits'
    wcs_hdu.writeto(output_filename, overwrite=True)

    #Delete temporary files. 
    os.remove(temp_filename)
    os.remove(astrometry_image_path)
    pdb.set_trace()