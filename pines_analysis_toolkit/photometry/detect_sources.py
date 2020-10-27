
import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
import pdb
from photutils import DAOStarFinder, aperture_photometry, CircularAperture
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
import pickle
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from pines_analysis_toolkit.data.bpm_chooser import bpm_chooser
from astropy.visualization import ImageNormalize, ZScaleInterval, SquaredStretch, SqrtStretch, LinearStretch
from pines_analysis_toolkit.data.bg_2d import bg_2d

'''Authors:
        Patrick Tamburo, Boston University, June 2020
	Purpose:
        Finds sources in reduced Mimir image. 
	Inputs:
        image_path (pathlib.Path): the path to the reduced image you want to detect sources in.
        source_dir (pathlib.Path): path to the source directory for this object for saving source detection image. 
        fhwm (float): the fwhm in pixels for the source detection.
        thresh (float): the number of sigma above background for source detection.
        plot (bool, optional): whether or not to plot the detected sources. Will also save plot. 
    Outputs:
        sources (astropy table): table of sources found in the image, containing centroids and quick photometry estimate.
	TODO:
        Choose fwhm for source detection based on seeing measured in log. 
        Make bpm a fits image?
'''

def detect_sources(image_path, source_dir, fwhm=5.0, thresh=2.0, plot=False):
    ap_rad = 5 #Radius of aperture in pixels for doing quick photometry on detected sources.
    edge_tolerance = 25 #Number of pixels from the edge where detected sources are cut from. We don't want these because they can shift off the detector.
    
    #Read in the image. 
    image = fits.open(image_path)[0].data
    header = fits.open(image_path)[0].header
    
    #Interpolate nans if any in image. 
    kernel = Gaussian2DKernel(x_stddev=1)
    image = interpolate_replace_nans(image, kernel)
    
    #Do a quick background model subtraction (makes source detection easier).
    image = bg_2d(image, box_size=64)

    #Get a bad pixel mask to help in source detection.
    pines_path = pines_dir_check()
    bpm_path = pines_path/('Calibrations/Bad Pixel Masks')            
    bpm, bad_pixel_mask_name = bpm_chooser(bpm_path, header)

    #Get the sigma_clipped_stats for the image.
    avg, med, std = sigma_clipped_stats(image)

    norm = ImageNormalize(data=image, interval=ZScaleInterval(), stretch=LinearStretch())

    if plot:
        title = image_path.name
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,9))
        ax.set_aspect('equal')
        im = ax.imshow(image, origin='lower', norm=norm, cmap='Greys_r')
        cax = fig.add_axes([0.87, 0.15, 0.035, 0.7])
        fig.colorbar(im, cax=cax, orientation='vertical', label='ADU')
        ax.set_title(title)

    print('')
    print('Finding sources in {}.'.format(image_path.name))
    print('....')
    print('')

    #Detect sources using DAOStarFinder.
    daofind = DAOStarFinder(fwhm=fwhm, threshold=thresh*std)  
    initial_sources = daofind(image - med, mask=bpm)

    #Do a cut based on source sharpness to get rid of some false detections.
    initial_sources.sort('sharpness')        
    bad_sharpness_locs = np.where(initial_sources['sharpness'] < 0.3)[0]
    initial_sources.remove_rows(bad_sharpness_locs)

    #Cut sources that are found within edge_tolerance pix of the edges.
    bad_x = np.where((initial_sources['xcentroid'] < edge_tolerance) | (initial_sources['xcentroid'] > 1023-edge_tolerance))[0]
    initial_sources.remove_rows(bad_x)
    bad_y = np.where((initial_sources['ycentroid'] < edge_tolerance) | (initial_sources['ycentroid'] > 1023-edge_tolerance))[0]
    initial_sources.remove_rows(bad_y)

    #Cut sources near y = 512, these are frequently bad. 
    bad_512 = np.where((initial_sources['ycentroid'] > 508) & (initial_sources['ycentroid'] < 516))
    initial_sources.remove_rows(bad_512)

    #Do quick photometry on the remaining sources. 
    positions = [(initial_sources['xcentroid'][i], initial_sources['ycentroid'][i]) for i in range(len(initial_sources))]
    apertures = CircularAperture(positions, r=ap_rad)
    phot_table = aperture_photometry(image-med, apertures, mask=bpm)

    #Cut based on brightness.
    phot_table.sort('aperture_sum')
    cutoff = 1*std*np.pi*ap_rad**2
    bad_source_locs = np.where(phot_table['aperture_sum'] < cutoff)
    phot_table.remove_rows(bad_source_locs)
    initial_sources.remove_rows(bad_source_locs)

    if plot:
        #Plot detected sources. 
        ax.plot(phot_table['xcenter'],phot_table['ycenter'],'rx')
        plt.tight_layout()
        plt.savefig(source_dir/('source_detection.png'))

    print('Found {} sources.'.format(len(phot_table)))
    sources = phot_table[::-1].to_pandas() #Resort remaining sources so that the brightest are listed firsts. 
    return sources
    