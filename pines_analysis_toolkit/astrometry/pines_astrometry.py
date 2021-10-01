import pines_analysis_toolkit as pat 
from glob import glob 
from natsort import natsorted 
import numpy as np
from astropy.convolution import interpolate_replace_nans, Gaussian2DKernel
from astropy.io import fits 
import os 
import time
from pathlib import Path

def pines_astrometry(target, api_key, download_data=False):
    """Uploads reduced images for a target to astrometry.net, downloads solution image, and updates the image header with the astrometry.net wcs.

    :param target: long name of the target
    :type target: str
    :param api_key: he api key of your account on astrometry.net
    :type api_key: str
    :param download_data: whether or not to first download reduced data for the target from the PINES server, defaults to False
    :type download_data: bool, optional
    """
    pines_path = pat.utils.pines_dir_check()
    short_name = pat.utils.short_name_creator(target)

    if download_data:
        sftp = pat.utils.pines_login()
        pat.data.get_reduced_science_files(sftp, target)

    reduced_path = pines_path/('Objects/'+short_name+'/reduced/')
    reduced_files = np.array(natsorted(glob(str(reduced_path)+'/*red.fits')))

    kernel = Gaussian2DKernel(x_stddev=0.25)

    times = []
    for i in range(0, len(reduced_files)):
        t1 = time.time()
        filename = Path(reduced_files[i])
        print(short_name+', '+filename.name+', image '+str(i+1)+' of '+str(len(reduced_files)))

        #If the header does not have a HISTORY keyword (which is added by astrometry.net), process it. 
        header = fits.open(filename)[0].header

        if 'HISTORY' not in header:
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
                print('Astrometry solution failed!')
                #Add HISTORY keyword to header so it doesn't get processed in the future. 
                header['HISTORY'] = 'Astrometric solution failed.'
                fits.writeto(filename, original_image, header, overwrite=True)

            #Delete temporary files. 
            os.remove(temp_filename)
            os.remove(astrometry_image_path)
            time.sleep(1)
            t2 = time.time()
            times.append(t2-t1)
            print('Average completion time: {:3.1f} s.'.format(np.mean(times)))
            print('')

        #If the header DOES have a HISTORY keyword, skip it, it has already been processed. 
        else:
            print('Astrometric processing already complete for {}, skipping.'.format(filename.name))
            print('')