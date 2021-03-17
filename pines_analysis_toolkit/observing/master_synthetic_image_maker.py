import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter
import pdb
import scipy.optimize as opt
from scipy import signal
from astropy.modeling import models, fitting
from astropy.io import fits
from photutils import DAOStarFinder, aperture_photometry, CircularAperture, Background2D, MedianBackground
from astropy.stats import sigma_clipped_stats, SigmaClip
import pickle 
from pathlib import Path
import os 
import time
import datetime
import shutil
import sys
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.utils.quick_plot import quick_plot as qp 

def master_synthetic_image_creator(target, image_name, seeing=2.5, sigma_above_bg=5):
    '''Authors:
            Patrick Tamburo, Boston University, February 2021
        Purpose:
            Creates a master synthetic image for a PINES target by detecting sources in a reduced image of the field.
        Inputs:
            target (str): The target's full 2MASS name.
            image_name (str): The name of the reduced file (e.g., 20210204.420_red.fits).
            seeing (float): The FWHM seeing of the image in arcsec. By default, 2.5". 
            sigma_above_bg (float): The sigma above background used to rule in sources in daostarfinder. By default, 5.
        Outputs:
            Writes out a master synthetic image to PINES_analysis_toolkit/Calibrations/Master Synthetic Images/.
        TODO:
            Write the filename used to create the master image to the master image header. 
    '''


    def mimir_source_finder(image_path,sigma_above_bg,fwhm, exclude_lower_left=False):
        """Find sources in Mimir images."""
        
        np.seterr(all='ignore') #Ignore invalids (i.e. divide by zeros)

        #Find stars in the master image.
        avg, med, stddev = sigma_clipped_stats(image,sigma=3.0,maxiters=3) #Previously maxiters = 5!   
        daofind = DAOStarFinder(fwhm=fwhm, threshold=sigma_above_bg*stddev,sky=med,ratio=0.8)
        new_sources = daofind(image)
        x_centroids = new_sources['xcentroid']
        y_centroids = new_sources['ycentroid']
        sharpness = new_sources['sharpness']
        fluxes = new_sources['flux']
        peaks = new_sources['peak']

        #Cut sources that are found within 20 pix of the edges.
        use_x = np.where((x_centroids > 20) & (x_centroids < 1004))[0]
        x_centroids = x_centroids[use_x]
        y_centroids = y_centroids[use_x]
        sharpness = sharpness[use_x]
        fluxes = fluxes[use_x]
        peaks = peaks[use_x]
        use_y = np.where((y_centroids  > 20) & (y_centroids  < 1004))[0]
        x_centroids  = x_centroids [use_y]
        y_centroids  = y_centroids [use_y]
        sharpness = sharpness[use_y]
        fluxes = fluxes[use_y]
        peaks = peaks[use_y]

        #Also cut using sharpness, this seems to eliminate a lot of false detections.
        use_sharp = np.where(sharpness > 0.5)[0]
        x_centroids  = x_centroids [use_sharp]
        y_centroids  = y_centroids [use_sharp]
        sharpness = sharpness[use_sharp]
        fluxes = fluxes[use_sharp]
        peaks = peaks[use_sharp]

        if exclude_lower_left:
            #Cut sources in the lower left, if bars are present.
            use_ll =  np.where((x_centroids > 512) | (y_centroids > 512))
            x_centroids  = x_centroids [use_ll]
            y_centroids  = y_centroids [use_ll]
            sharpness = sharpness[use_ll]
            fluxes = fluxes[use_ll]
            peaks = peaks[use_ll]
        
        #Cut targets whose y centroids are near y = 512. These are usually bad.
        use_512 = np.where(np.logical_or((y_centroids < 510),(y_centroids > 514)))[0]
        x_centroids  = x_centroids [use_512]
        y_centroids  = y_centroids [use_512]
        sharpness = sharpness[use_512]
        fluxes = fluxes[use_512]
        peaks = peaks[use_512]

        #Cut sources with negative/saturated peaks
        use_peaks = np.where((peaks > 30) & (peaks < 3000))[0]
        x_centroids  = x_centroids [use_peaks]
        y_centroids  = y_centroids [use_peaks]
        sharpness = sharpness[use_peaks]
        fluxes = fluxes[use_peaks]
        peaks = peaks[use_peaks]

        #Do quick photometry on the remaining sources. 
        positions = [(x_centroids[i], y_centroids[i]) for i in range(len(x_centroids))]
        apertures = CircularAperture(positions, r=4)
        phot_table = aperture_photometry(image-med, apertures)

        #Cut based on brightness.
        phot_table.sort('aperture_sum')
        cutoff = 1*std*np.pi*4**2
        bad_source_locs = np.where(phot_table['aperture_sum'] < cutoff)
        phot_table.remove_rows(bad_source_locs)
        
        x_centroids = phot_table['xcenter'].value
        y_centroids = phot_table['ycenter'].value

        pdb.set_trace()
        return(x_centroids,y_centroids)

    def synthetic_image_maker(x_centroids,y_centroids,fwhm):
        #Construct synthetic images from centroid/flux data.
        synthetic_image = np.zeros((1024,1024))
        sigma = fwhm/2.355
        for i in range(len(x_centroids)):
            #Cut out little boxes around each source and add in Gaussian representations. This saves time. 
            int_centroid_x = int(np.round(x_centroids[i]))
            int_centroid_y = int(np.round(y_centroids[i]))
            y_cut, x_cut = np.mgrid[int_centroid_y-10:int_centroid_y+10,int_centroid_x-10:int_centroid_x+10]
            dist = np.sqrt((x_cut-x_centroids[i])**2+(y_cut-y_centroids[i])**2)
            synthetic_image[y_cut,x_cut] += np.exp(-((dist)**2/(2*sigma**2)+((dist)**2/(2*sigma**2))))
        return(synthetic_image)
    
    target = target.replace(' ','')
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    master_synthetic_path = pines_path/('Calibrations/Master Synthetic Images/'+target+'_master_synthetic.fits')
    image_path = pines_path/('Objects/'+short_name+'/reduced/'+image_name)
    plt.ion()

    target = target.replace(' ','')
    seeing = float(seeing)
    daostarfinder_fwhm = seeing*2.355/0.579


    #Open the image and calibration files. 
    header = fits.open(image_path)[0].header
    image = fits.open(image_path)[0].data

    #Interpolate over bad pixels
    kernel = Gaussian2DKernel(x_stddev=1)
    image = interpolate_replace_nans(image, kernel)

    #Do a simple 2d background model. 
    box_size=32
    sigma_clip = SigmaClip(sigma=3.)
    bkg_estimator = MedianBackground()
    bkg = Background2D(image, (box_size, box_size), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    image = image - bkg.background

    avg,med,std = sigma_clipped_stats(image)

    #Find sources in the image. 
    (x_centroids,y_centroids) = mimir_source_finder(image, sigma_above_bg=sigma_above_bg, fwhm=daostarfinder_fwhm)

    #Plot the field with detected sources. 
    qp(image)
    plt.plot(x_centroids,y_centroids,'rx')
    for i in range(len(x_centroids)):
        plt.text(x_centroids[i]+8,y_centroids[i]+8,str(i),color='r', fontsize=14)
    plt.title('Inspect to make sure stars were found!\nO for magnification tool, R to reset view')
    plt.tight_layout()
    plt.show()

    print('')
    print('')
    print('')

    #Prompt the user to remove any false detections.
    ids = input('Enter ids of sources to be removed separated by commas (i.e., 4,18,22). If none to remove, hit enter. To break, ctrl + D. ')
    if ids != '':
        ids_to_eliminate = [int(i) for i in ids.split(',')]
        ids = [int(i) for i in np.linspace(0,len(x_centroids)-1,len(x_centroids))]
        ids_to_keep = []
        for i in range(len(ids)):
            if ids[i] not in ids_to_eliminate:
                ids_to_keep.append(ids[i])
    else:
        ids_to_keep = [int(i) for i in np.linspace(0,len(x_centroids)-1,len(x_centroids))]
    plt.clf()
    plt.imshow(image,origin='lower',vmin=med,vmax=med+5*std)
    plt.plot(x_centroids[ids_to_keep],y_centroids[ids_to_keep],'rx')
    for i in range(len(x_centroids[ids_to_keep])):
        plt.text(x_centroids[ids_to_keep][i]+8,y_centroids[ids_to_keep][i]+8,str(i),color='r')

    #Create the synthetic image using the accepted sources. 
    synthetic_image = synthetic_image_maker(x_centroids[ids_to_keep],y_centroids[ids_to_keep],8)
    plt.figure(figsize=(9,7))
    plt.imshow(synthetic_image,origin='lower')
    plt.title('Synthetic image')
    plt.show()

    pdb.set_trace()
    print('')
    print('')
    print('')
    #Now write to a master synthetic image.fits file.
    hdu = fits.PrimaryHDU(synthetic_image)
    hdu.writeto(master_synthetic_path, overwrite=True)
    print('Writing master synthetic image to {}/'.format(master_synthetic_path))


if __name__ == '__main__':
    target = '2MASS J10292165+1626526'
    image_name = '20210303.197_red.fits'
    seeing = 2.3
    master_synthetic_image_creator(target, image_name, seeing=seeing, sigma_above_bg=10)