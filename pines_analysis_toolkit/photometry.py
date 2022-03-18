import pines_analysis_toolkit as pat
from pines_analysis_toolkit.analysis import raw_flux_plot
from pines_analysis_toolkit.utils import pines_dir_check, short_name_creator, pines_log_reader, quick_plot as qp, jd_utc_to_bjd_tdb, get_source_names, jd_utc_to_bjd_tdb
from pines_analysis_toolkit.data import master_dark_stddev_chooser, bpm_chooser, bg_2d
from pines_analysis_toolkit.analysis import night_splitter, block_splitter
from pines_analysis_toolkit.observing import seeing_measurer

from astropy.io import fits
from astropy.stats import sigma_clipped_stats, gaussian_sigma_to_fwhm
from astropy.table import Table, QTable
import astropy.units as u
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans, convolve, Gaussian1DKernel
from astropy.visualization import ZScaleInterval, ImageNormalize, SquaredStretch, LinearStretch, SqrtStretch, MinMaxInterval, simple_norm
from astropy.modeling import models, fitting
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.nddata import NDData
from astropy.utils.exceptions import AstropyUserWarning
from astropy.coordinates import SkyCoord, Angle

from photutils import CircularAperture, CircularAnnulus, aperture_photometry, make_source_mask, centroid_2dg, centroid_1dg, centroid_com, EPSFBuilder, aperture_photometry, EPSFFitter
from photutils.psf import BasicPSFPhotometry, IntegratedGaussianPRF, DAOGroup, extract_stars, IterativelySubtractedPSFPhotometry
from photutils.utils import calc_total_error
from photutils.detection import IRAFStarFinder, DAOStarFinder
from photutils.centroids import centroid_sources
from photutils.background import MMMBackground

from scipy.stats import sigmaclip, mode as scipy_mode
from scipy import optimize
from scipy import interpolate
from scipy.optimize import least_squares


import numpy as np
from natsort import natsorted
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import math
import shutil
import warnings
import os
import time
from glob import glob
import julian

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.pyplot import Circle
from progressbar import ProgressBar
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
import matplotlib
matplotlib.use('TkAgg')


def hmsm_to_days(hour=0, min=0, sec=0, micro=0):
    """
    Convert hours, minutes, seconds, and microseconds to fractional days.

    """
    days = sec + (micro / 1.e6)
    days = min + (days / 60.)
    days = hour + (days / 60.)
    return days / 24.


def date_to_jd(year, month, day):
    """
    Convert a date to Julian Day.

    Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet', 
        4th ed., Duffet-Smith and Zwart, 2011.

    """
    if month == 1 or month == 2:
        yearp = year - 1
        monthp = month + 12
    else:
        yearp = year
        monthp = month

    # this checks where we are in relation to October 15, 1582, the beginning
    # of the Gregorian calendar.
    if ((year < 1582) or
        (year == 1582 and month < 10) or
            (year == 1582 and month == 10 and day < 15)):
        # before start of Gregorian calendar
        B = 0
    else:
        # after start of Gregorian calendar
        A = math.trunc(yearp / 100.)
        B = 2 - A + math.trunc(A / 4.)

    if yearp < 0:
        C = math.trunc((365.25 * yearp) - 0.75)
    else:
        C = math.trunc(365.25 * yearp)

    D = math.trunc(30.6001 * (monthp + 1))

    jd = B + C + D + day + 1720994.5

    return jd


def iraf_style_photometry(phot_apertures, bg_apertures, data, dark_std_data, header, seeing, bg_method='mean', epadu=1.0, gain=8.21, non_linear_threshold=4000):
    """Computes photometry with PhotUtils apertures, with IRAF formulae
    Parameters
    ----------
    phot_apertures : photutils PixelAperture object (or subclass)
        The PhotUtils apertures object to compute the photometry.
        i.e. the object returned via CircularAperture.
    bg_apertures : photutils PixelAperture object (or subclass)
        The phoutils aperture object to measure the background in.
        i.e. the object returned via CircularAnnulus.
    data : array
        The data for the image to be measured.
    bg_method: {'mean', 'median', 'mode'}, optional
        The statistic used to calculate the background.
        All measurements are sigma clipped.
        NOTE: From DAOPHOT, mode = 3 * median - 2 * mean.
    epadu: float, optional
        Gain in electrons per adu (only use if image units aren't e-).
    Returns
    -------
    final_tbl : astropy.table.Table
        An astropy Table with the colums X, Y, flux, flux_error, mag,
        and mag_err measurements for each of the sources.
    """
    exptime = header['EXPTIME']

    if bg_method not in ['mean', 'median', 'mode']:
        raise ValueError(
            'Invalid background method, choose either mean, median, or mode')

    # Create a list to hold the flux for each source.
    aperture_sum = []
    interpolation_flags = np.zeros(len(phot_apertures.positions), dtype='bool')
    
    #TESTING, REMOVE LATER!!!
    target_pixel = (421, 864)
    save_val = data[target_pixel[1], target_pixel[0]]

    
    # #NEW!! Joint PSF fitting approach
    # #Create an array of all source cutouts for this frame
    # cutout_array = cutout_array_creator(data, phot_apertures.positions[:,0], phot_apertures.positions[:,1], cutout_w=11)
    # #Pass the cutout array to the joint psf fitter
    # joint_fit_params, xs, ys = joint_psf_fitter(cutout_array)

    for i in range(len(phot_apertures.positions)):
        pos = phot_apertures.positions[i]
        # Cutout around the source position
        cutout_w = 30
        x_pos = pos[0]
        y_pos = pos[1]
        cutout = data[int((y_pos-cutout_w)):int(y_pos+cutout_w)+1, int(x_pos-cutout_w):int(x_pos+cutout_w)+1]



        x_cutout = x_pos - np.floor(x_pos - cutout_w)
        y_cutout = y_pos - np.floor(y_pos - cutout_w)
        ap = CircularAperture((x_cutout, y_cutout), r=phot_apertures.r)

        if i == 3:
            
            an = CircularAnnulus((x_cutout, y_cutout), r_in=12, r_out=30)
            mask = an.to_mask(method='center')
            mask_pix = mask.multiply(cutout)
            an_pix = mask_pix[mask_pix != 0]
            an_pix = an_pix[~np.isnan(an_pix)]
            v, l, h = sigmaclip(an_pix)
            sky = np.mean(v)
            star = save_val - sky
            noise_estimate = np.sqrt(star + sky + (22.65/exptime) + (26.6**2))


            #Change the target pixel's value using random noise estimate.  
            data[target_pixel[1], target_pixel[0]] = np.random.normal(data[target_pixel[1], target_pixel[0]], noise_estimate)

            #Rebuild the cutout. 
            cutout = data[int((y_pos-cutout_w)):int(y_pos+cutout_w)+1, int(x_pos-cutout_w):int(x_pos+cutout_w)+1]


        # Cut out the pixels JUST inside the aperture, and check if there are NaNs there. If so, interpolate over NaNs in the cutout.
        ap_mask = ap.to_mask(method='exact')
        ap_cut = ap_mask.cutout(cutout)

        try:
            bad_sum = np.sum(np.isnan(ap_cut))
        except:
            breakpoint
        
        # #NEW!! Joint PSF fitting approach
        # if bad_sum > 0:
        #     bads = np.where(np.isnan(cutout))
            
        #     joint_psf_cutout_offset = int(cutout_w - np.shape(cutout_array)[1]/2)
        #     bad_pixels_in_joint_psf =  np.where((bads[0]-joint_psf_cutout_offset >= 0) & (bads[0]-joint_psf_cutout_offset < np.shape(cutout_array)[1]) & (bads[1]-joint_psf_cutout_offset >= 0) & (bads[1]-joint_psf_cutout_offset < np.shape(cutout_array)[1]))[0]
            
        #     bad_dists = np.sqrt((bads[0] - y_cutout)** 2 + (bads[1] - x_cutout)**2)
        #     if np.sum(bad_dists < phot_apertures.r+1):
                
        #         #Generate a model using the joint fit
        #         model = twoD_Gaussian(xs[i,:,:], ys[i,:,:], joint_fit_params['sigma_x'], joint_fit_params['sigma_y'], joint_fit_params['rho'], joint_fit_params['xcens'][i], joint_fit_params['ycens'][i], joint_fit_params['amplitudes'][i], joint_fit_params['offsets'][i])
                
        #         #Replace the bad pixels that are within the joint psf cutout array with their values from the model
        #         #cutout[bads] = model[bads]
        #         if i == 3:
        #             print(save_val, model[bads[0][bad_pixels_in_joint_psf]-joint_psf_cutout_offset, bads[1][bad_pixels_in_joint_psf]-joint_psf_cutout_offset])
        #         cutout[bads[0][bad_pixels_in_joint_psf], bads[1][bad_pixels_in_joint_psf]] = model[bads[0][bad_pixels_in_joint_psf]-joint_psf_cutout_offset, bads[1][bad_pixels_in_joint_psf]-joint_psf_cutout_offset]
            
       #OLD!! 2D Gaussian approach
        if bad_sum > 0:
           bads = np.where(np.isnan(cutout))
           bad_dists = np.sqrt((bads[0] - y_cutout)** 2 + (bads[1] - x_cutout)**2)

           # Check if any bad pixels fall within the aperture. If so, set the interpolation flag to True for this source.
           if np.sum(bad_dists < phot_apertures.r+1):

               # 2D gaussian fitting approach
               # Set up a 2D Gaussian model to interpolate the bad pixel values in the cutout.
               model_init = models.Const2D(amplitude=np.nanmedian(cutout))+models.Gaussian2D(
                   amplitude=np.nanmax(cutout), x_mean=x_cutout, y_mean=y_cutout, x_stddev=seeing, y_stddev=seeing)
               # 2D grids of x and y coordinates
               xx, yy = np.indices(cutout.shape)
               # Find locations where the cutout has *good* values (i.e. not NaNs).
               mask = ~np.isnan(cutout)
               x = xx[mask]  # Only use coordinates at these good values.
               y = yy[mask]
               # Make a 1D cutout using only the good values.
               cutout_1d = cutout[mask]
               fitter = fitting.LevMarLSQFitter()
               # Fit the model to the 1d cutout.
               model_fit = fitter(model_init, x, y, cutout_1d)

               # Interpolate the pixels in the cutout using the 2D Gaussian fit.
               cutout[~mask] = model_fit(xx, yy)[~mask]

               interpolation_flags[i] = True

               # #Gaussian convolution approach.
               # cutout = interpolate_replace_nans(cutout, kernel=Gaussian2DKernel(x_stddev=0.5))
        
        phot_source = aperture_photometry(cutout, ap)
        # if np.isnan(phot_source['aperture_sum'][0]):
        #     breakpoint()

        aperture_sum.append(phot_source['aperture_sum'][0])

    # Add positions/fluxes to a table
    xcenter = phot_apertures.positions[:, 0]*u.pix
    ycenter = phot_apertures.positions[:, 1]*u.pix
    phot = QTable([xcenter, ycenter, aperture_sum], names=('xcenter', 'ycenter', 'aperture_sum'))

    # Now measure the background around each source.
    # Make a mask to block out any sources that might show up in the annuli and bias them.
    mask = make_source_mask(data, nsigma=3, npixels=5, dilate_size=7)
    # Pass the data with sources masked out to the bg calculator.
    bg_phot = aperture_stats_tbl(~mask*data, bg_apertures, sigma_clip=True)
    ap_area = phot_apertures.area

    bg_method_name = 'aperture_{}'.format(bg_method)
    background = bg_phot[bg_method_name]

    flux = phot['aperture_sum'] - background * ap_area

    # Need to use variance of the sources for Poisson noise term in error computation.
    flux_error = compute_phot_error(
        flux, bg_phot, bg_method, ap_area, exptime, dark_std_data, phot_apertures, epadu)

    mag = -2.5 * np.log10(flux)
    mag_err = 1.0857 * flux_error / flux

    # Make the final table
    X, Y = phot_apertures.positions.T
    stacked = np.stack([X, Y, flux, flux_error, mag, mag_err,
                       background, interpolation_flags], axis=1)
    names = ['X', 'Y', 'flux', 'flux_error', 'mag',
             'mag_error', 'background', 'interpolation_flag']

    final_tbl = Table(data=stacked, names=names)

    # Check for nans
    if sum(np.isnan(final_tbl['flux'])) > 0:
        bad_locs = np.where(np.isnan(final_tbl['flux']))[0]
        # breakpoint()

    return final_tbl


def compute_phot_error(flux, bg_phot, bg_method, ap_area, exptime, dark_std_data, phot_apertures, epadu=1.0):
    """Computes the flux errors using the DAOPHOT style computation.
        Includes photon noise from the source, background terms, dark current, and read noise."""

    # See eqn. 1 from Broeg et al. (2005): https://ui.adsabs.harvard.edu/abs/2005AN....326..134B/abstract
    flux_variance_term = flux / epadu  # This is just the flux in photons.
    bg_variance_term_1 = ap_area * (bg_phot['aperture_std'])**2
    bg_variance_term_2 = (
        ap_area * bg_phot['aperture_std'])**2 / bg_phot['aperture_area']

    # Measure combined read noise + dark current using the dark standard deviation image.
    # TODO: Check that this is correct.
    dark_rn_term = np.zeros(len(flux_variance_term))
    for i in range(len(phot_apertures)):
        ap = phot_apertures[i]
        dark_rn_ap = ap.to_mask().multiply(dark_std_data)
        dark_rn_ap = dark_rn_ap[dark_rn_ap != 0]
        dark_rn_term[i] = ap_area * (np.median(dark_rn_ap)**2)

    variance = flux_variance_term + bg_variance_term_1 + \
        bg_variance_term_2 + dark_rn_term
    flux_error = variance ** .5
    return flux_error


def aperture_stats_tbl(data, apertures, method='center', sigma_clip=True):
    """Computes mean/median/mode/std in Photutils apertures.
    Compute statistics for custom local background methods.
    This is primarily intended for estimating backgrounds
    via annulus apertures.  The intent is that this falls easily
    into other code to provide background measurements.
    Parameters
    ----------
    data : array
        The data for the image to be measured.
    apertures : photutils PixelAperture object (or subclass)
        The phoutils aperture object to measure the stats in.
        i.e. the object returned via CirularAperture,
        CircularAnnulus, or RectangularAperture etc.
    method: str
        The method by which to handle the pixel overlap.
        Defaults to computing the exact area.
        NOTE: Currently, this will actually fully include a
        pixel where the aperture has ANY overlap, as a median
        is also being performed.  If the method is set to 'center'
        the pixels will only be included if the pixel's center
        falls within the aperture.
    sigma_clip: bool
        Flag to activate sigma clipping of background pixels
    Returns
    -------
    stats_tbl : astropy.table.Table
        An astropy Table with the colums X, Y, aperture_mean,
        aperture_median, aperture_mode, aperture_std, aperture_area
        and a row for each of the positions of the apertures.
    """

    # Get the masks that will be used to identify our desired pixels.
    masks = apertures.to_mask(method=method)

    # Compute the stats of pixels within the masks
    aperture_stats = [calc_aperture_mmm(data, mask, sigma_clip) for mask in masks]

    aperture_stats = np.array(aperture_stats)

    # Place the array of the x y positions alongside the stats
    stacked = np.hstack([apertures.positions, aperture_stats])

    # Name the columns
    names = ['X', 'Y', 'aperture_mean', 'aperture_median',
             'aperture_mode', 'aperture_std', 'aperture_area']

    # Make the table
    stats_tbl = Table(data=stacked, names=names)

    return stats_tbl


def calc_aperture_mmm(data, mask, sigma_clip):
    """Helper function to actually calculate the stats for pixels falling within some Photutils aperture mask on some array of data."""
    #cutout = mask.cutout(data, fill_value=np.nan)
    cutout = mask.multiply(data)
    if cutout is None:
        return (np.nan, np.nan, np.nan, np.nan, np.nan)
    else:
        values = cutout[cutout != 0]  # Unwrap the annulus into a 1D array.
        values = values[~np.isnan(values)]  # Ignore any NaNs in the annulus.
        if sigma_clip:
            # Sigma clip the 1D annulus values.
            values, clow, chigh = sigmaclip(values, low=3, high=3)
        mean = np.nanmean(values)
        median = np.nanmedian(values)
        std = np.nanstd(values)
        mode = 3 * median - 2 * mean
        actual_area = (~np.isnan(values)).sum()

        return (mean, median, mode, std, actual_area)


def fixed_aper_phot(short_name, centroided_sources, ap_radii, an_in=12., an_out=30., gain=8.21, qe=0.9, force_output_path=''):
    '''Authors:
                Patrick Tamburo, Boston University, June 2020
        Purpose:
        Performs *fixed* aperture photometry on a set of reduced images given dataframe of source positions.
        The iraf_style_photometry, compute_phot_error, perture_stats_tbl, and calc_aperture_mmm routines are from Varun Bajaj on github:
            https://github.com/spacetelescope/wfc3_photometry/blob/master/photometry_tools/photometry_with_errors.py. 
        Inputs:
        target (str): The target's full 2MASS name.
        sources (pandas dataframe): List of source names, x and y positions in every image. 
        ap_radii (list of floats): List of aperture radii in pixels for which aperture photometry wil be performed. 
        an_in (float, optional): The inner radius of the annulus used to estimate background, in pixels. 
        an_out (float, optional): The outer radius of the annulus used to estimate background, in pixels. 
        gain (float, optional): The gain of the detector in e-/ADU.
        qe (float, optional): The quantum efficiency of the detector.
        force_output_path (path): if you want to manually set an output directory, specify the top-level here (i.e. the directory containing Logs, Objects, etc.)
    Outputs:
        Saves aperture photometry csv to PINES_analysis_toolkit/Objects/short_name/aper_phot/ for each aperture.
        TODO:
    '''

    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()

    # Remove any leading/trailing spaces in the column names.
    centroided_sources.columns = centroided_sources.columns.str.lstrip()
    centroided_sources.columns = centroided_sources.columns.str.rstrip()

    # Get list of reduced files for target.
    reduced_path = pines_path/('Objects/'+short_name+'/reduced')
    reduced_filenames = natsorted(
        [x.name for x in reduced_path.glob('*red.fits')])
    reduced_files = np.array([reduced_path/i for i in reduced_filenames])

    source_names = get_source_names(centroided_sources)

    
     # Declare a new dataframe to hold the information for all targets for this aperture.
    columns = ['Filename', 'Time UT', 'Time JD UTC', 'Time BJD TDB',
                'Night Number', 'Block Number', 'Airmass', 'Seeing']
    for i in range(0, len(source_names)):
        columns.append(source_names[i]+' Flux')
        columns.append(source_names[i]+' Flux Error')
        columns.append(source_names[i]+' Background')
        columns.append(source_names[i]+' Interpolation Flag')

    ap_dfs = []
    for i in range(len(ap_radii)):
        ap_dfs.append(pd.DataFrame(index=range(len(reduced_files)), columns=columns))
    
    print('Doing fixed aperture photometry with radii = {} pixels.'.format(ap_radii))
    pbar = ProgressBar()
    for j in pbar(range(len(reduced_files))):
        hdu = fits.open(reduced_files[j])[0]
        data = hdu.data
        header = hdu.header

        # Get the source positions in this image.
        positions = []
        for i in range(len(source_names)):
            positions.append((float(centroided_sources[source_names[i]+' Image X'][j]), float(centroided_sources[source_names[i]+' Image Y'][j])))

        # Read in some supporting information.
        log_path = pines_path / ('Logs/'+reduced_files[j].name.split('.')[0]+'_log.txt')
        log = pines_log_reader(log_path)
        log_ind = np.where(log['Filename'] == reduced_files[j].name.split('_')[0]+'.fits')[0][0]
        
        # Get the closest date master_dark_stddev image for this exposure time.
        # We'll use this to measure read noise and dark current.
        master_dark_stddev = master_dark_stddev_chooser(pines_path/('Calibrations/Darks/Master Darks Stddev/'), header)

        date = pat.utils.mimir_date_reader(header)
        jd = julian.to_jd(date)
        filename = reduced_files[j].name
        time_ut = header['DATE-OBS']
        time_bjd = jd_utc_to_bjd_tdb(jd, header['TELRA'], header['TELDEC'])
        night_number = centroided_sources['Night Number'][j]
        block_number = centroided_sources['Block Number'][j]
        airmass = header['AIRMASS']
        seeing = log['X seeing'][log_ind]
        for i in range(len(ap_radii)):
            ap_df = ap_dfs[i] 
            ap_df['Filename'][j] = filename
            ap_df['Time UT'][j] = time_ut
            ap_df['Time JD UTC'][j] = jd
            # Using the telescope ra and dec should be accurate enough for our purposes
            ap_df['Time BJD TDB'][j] = time_bjd
            ap_df['Night Number'][j] = night_number
            ap_df['Block Number'][j] = block_number
            ap_df['Airmass'][j] = airmass
            ap_df['Seeing'][j] = seeing

        # If the shift quality has been flagged, skip this image.
        if log['Shift quality flag'].iloc[log_ind] == 1:
            continue

        for i in range(len(ap_radii)):
            ap_df = ap_dfs[i]

            # Create an aperture centered on this position with radius = ap.
            try:
                apertures = CircularAperture(positions, r=ap_radii[i])
            except:
                breakpoint()

            # Create an annulus centered on this position.
            annuli = CircularAnnulus(positions, r_in=an_in, r_out=an_out)

            photometry_tbl = iraf_style_photometry(apertures, annuli, data*gain, master_dark_stddev*gain, header, ap_df['Seeing'][j])

            for k in range(len(source_names)):
                ap_df[source_names[k]+ ' Flux'][j] = photometry_tbl['flux'][k]
                ap_df[source_names[k] + ' Flux Error'][j] = photometry_tbl['flux_error'][k]
                ap_df[source_names[k] + ' Background'][j] = photometry_tbl['background'][k]
                ap_df[source_names[k]+ ' Interpolation Flag'][j] = int(photometry_tbl['interpolation_flag'][k])
    
    for i in range(len(ap_radii)):
        ap_df = ap_dfs[i]
        ap = ap_radii[i] 
        output_filename = pines_path/('Objects/'+short_name+'/aper_phot/' + short_name+'_fixed_aper_phot_{:1.1f}_pix_radius.csv'.format(float(ap)))

        # Write output to file.
        print('Saving ap = {:1.1f} aperture photometry output to {}.'.format(
            ap, output_filename))
        print('')
        with open(output_filename, 'w') as f:
            for j in range(len(ap_df)):
                # Write in the header.
                if j == 0:
                    f.write('{:>21s}, {:>22s}, {:>17s}, {:>17s}, {:>17s}, {:>17s}, {:>7s}, {:>7s}, '.format(
                        'Filename', 'Time UT', 'Time JD UTC', 'Time BJD TDB', 'Night Number', 'Block Number', 'Airmass', 'Seeing'))
                    for i in range(len(source_names)):
                        if i != len(source_names) - 1:
                            f.write('{:>22s}, {:>28s}, {:>28s}, {:>34s}, '.format(
                                source_names[i]+' Flux', source_names[i]+' Flux Error', source_names[i]+' Background', source_names[i]+' Interpolation Flag'))
                        else:
                            f.write('{:>22s}, {:>28s}, {:>28s}, {:>34s}\n'.format(
                                source_names[i]+' Flux', source_names[i]+' Flux Error', source_names[i]+' Background', source_names[i]+' Interpolation Flag'))

                # Write in Filename, Time UT, Time JD, Airmass, Seeing values.
                format_string = '{:21s}, {:22s}, {:17.9f}, {:17.9f}, {:17d}, {:17d}, {:7.2f}, {:7.1f}, '
                # If the seeing value for this image is 'nan' (a string), convert it to a float.
                # TODO: Not sure why it's being read in as a string, fix that.
                if type(ap_df['Seeing'][j]) == str:
                    ap_df['Seeing'][j] = float(ap_df['Seeing'][j])

                # Do a try/except clause for writeout, in case it breaks in the future.
                try:
                    f.write(format_string.format(ap_df['Filename'][j], ap_df['Time UT'][j], ap_df['Time JD UTC'][j], ap_df['Time BJD TDB']
                            [j], ap_df['Night Number'][j], ap_df['Block Number'][j], ap_df['Airmass'][j], ap_df['Seeing'][j]))
                except:
                    print(
                        'Writeout failed! Inspect quantities you are trying to write out.')
                    breakpoint()

                # Write in Flux, Flux Error, and Background values for every source.
                for i in range(len(source_names)):
                    if i != len(source_names) - 1:
                        format_string = '{:22.5f}, {:28.5f}, {:28.5f}, {:34d}, '
                    else:
                        format_string = '{:22.5f}, {:28.5f}, {:28.5f}, {:34d}\n'
                    try:
                        f.write(format_string.format(ap_df[source_names[i]+' Flux'][j], ap_df[source_names[i]+' Flux Error']
                                [j], ap_df[source_names[i]+' Background'][j], ap_df[source_names[i]+' Interpolation Flag'][j]))
                    except:
                        if i != len(source_names) - 1:
                            format_string = '{:22.5f}, {:28.5f}, {:28.5f}, {:34f}, '
                        else:
                            format_string = '{:22.5f}, {:28.5f}, {:28.5f}, {:34f}\n'
                        f.write(format_string.format(ap_df[source_names[i]+' Flux'][j], ap_df[source_names[i]+' Flux Error']
                                [j], ap_df[source_names[i]+' Background'][j], ap_df[source_names[i]+' Interpolation Flag'][j]))

        raw_flux_plot(output_filename, mode='night')
        print('')
    return

def variable_aper_phot(short_name, centroided_sources, multiplicative_factors, an_in=12., an_out=30., plots=False, gain=8.21, qe=0.9, plate_scale=0.579, force_output_path='', smoothing_size=5, bin_mins=0.0, time_threshold=0.1):
    plt.ioff()

    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()

    # Remove any leading/trailing spaces in the column names.
    centroided_sources.columns = centroided_sources.columns.str.lstrip()
    centroided_sources.columns = centroided_sources.columns.str.rstrip()

    # Get list of reduced files for target.
    reduced_path = pines_path/('Objects/'+short_name+'/reduced')
    reduced_filenames = natsorted(
        [x.name for x in reduced_path.glob('*red.fits')])
    reduced_files = np.array([reduced_path/i for i in reduced_filenames])

    # Get source names.
    source_names = get_source_names(centroided_sources)

    # Get seeing from the logs. 
    # night_inds = night_splitter(centroided_sources['Time BJD TDB'])
    # log_path = pines_path/'Logs/'
    # log_times = []
    # log_seeing = []
    # for i in range(len(night_inds)):
    #     log = log_path/(centroided_sources['Filename'][night_inds[i][0]].split('.')[0]+'_log.txt')
    #     log_df = pines_log_reader(log)
    #     dates = np.array(log_df['Date'])
    #     for j in range(len(dates)):
    #         if len(dates[j].split(':')[-1].split('.')[0]) == 3:
    #             dates[j] = dates[j].split(':')[0] + ':' + dates[j].split(':')[1] + ':' + dates[j].split(':')[-1][1:]

    #         if dates[j].split(':')[-1] == '60.00':
    #             dates[j] = dates[j].split(
    #                 ':')[0]+':'+str(int(dates[j].split(':')[1])+1)+':00.00'
    #         # Keep a try/except clause here in case other unknown DATE-OBS formats pop up.
    #         if '.' not in dates[j]:
    #             date_fmt = '%Y-%m-%dT%H:%M:%S'
    #         else:
    #             date_fmt = '%Y-%m-%dT%H:%M:%S.%f'
    #         date = datetime.datetime.strptime(dates[j], date_fmt)
    #     jd_utc = julian.to_jd(date)

    #     breakpoint()
    #     log_times.extend()
    #     breakpoint()
    # seeing = np.array(centroided_sources['Seeing'], dtype='float')

    seeing_csv_path = pines_path/('Objects/'+short_name+'/sources/seeing.csv')
    if not os.path.exists(seeing_csv_path):
        print('Seeing .csv does not exist! Using log measurements...')
        #seeing_measurer(pines_path/('Objects/'+short_name))
        seeing = np.array(centroided_sources['Seeing'], dtype='float')
    else:
        seeing_df = pines_log_reader(seeing_csv_path)
        seeing = np.array(seeing_df['Seeing'], dtype='float')
    
    seeing[np.isnan(seeing)] = np.nanmedian(seeing)

    #Get block indices for setting up smoothed seeing. 
    block_inds = block_splitter(centroided_sources['Time BJD TDB'], bin_mins=bin_mins, time_threshold=time_threshold)
    #Get running average of seeing. 
    #Force the size to be odd. 
    if smoothing_size % 2 == 0:
        smoothing_size -= 1

#    plt.figure(figsize=(16,5))
    moving_avg_seeing = []
    for i in range(len(block_inds)):
        smoothed_block_seeing = np.convolve(seeing[block_inds[i]], np.ones(smoothing_size)/smoothing_size, mode='same')
        #Fix the edges. 
        smoothed_block_seeing[0:int((smoothing_size+1)/2)] = np.mean(smoothed_block_seeing[int((smoothing_size+1)/2):2*int((smoothing_size+1)/2)])
        smoothed_block_seeing[-int((smoothing_size+1)/2):] = np.mean(smoothed_block_seeing[-2*int((smoothing_size+1)/2):-int((smoothing_size+1)/2)])
        #Smooth with a Gaussian kernel
        smoothed_block_seeing = convolve(smoothed_block_seeing, Gaussian1DKernel(1.5), boundary='extend')
        moving_avg_seeing.extend(smoothed_block_seeing)
    
        # plt.plot(centroided_sources['Time BJD TDB'][block_inds[i]], seeing[block_inds[i]]/plate_scale, marker='.', alpha=0.5, color='tab:blue')
        # plt.plot(centroided_sources['Time BJD TDB'][block_inds[i]], smoothed_block_seeing/plate_scale, lw=2, color='tab:orange')
    #plt.tight_layout()
    moving_avg_seeing = np.array(moving_avg_seeing)

    # Loop over multiplicative factors
    for i in range(len(multiplicative_factors)):
        fact = multiplicative_factors[i]
        if np.max(moving_avg_seeing) * fact / plate_scale > an_in:
            an_in = int(np.max(moving_avg_seeing) * fact / plate_scale)+1
            an_out = int(np.max(moving_avg_seeing) * fact / plate_scale)+19

            print('WARNING: this will use apertures larger than inner annulus radius, extending it outward to {} pixels.'.format(an_in))

        print('Doing variable aperture photometry for {}, multiplicative seeing factor = {}, inner annulus radius = {} pix, outer annulus radius = {} pix.'.format(short_name, fact, an_in, an_out))

        # Declare a new dataframe to hold the information for all targets for this aperture.
        columns = ['Filename', 'Time UT', 'Time JD UTC', 'Time BJD TDB',
                   'Night Number', 'Block Number', 'Airmass', 'Seeing']
        for j in range(0, len(source_names)):
            columns.append(source_names[j]+' Flux')
            columns.append(source_names[j]+' Flux Error')
            columns.append(source_names[j]+' Background')
            columns.append(source_names[j]+' Interpolation Flag')

        var_df = pd.DataFrame(index=range(len(reduced_files)), columns=columns)
        output_filename = pines_path/('Objects/'+short_name+'/aper_phot/' + short_name+'_variable_aper_phot_{:1.2f}_seeing_factor.csv'.format(fact))

        # plt.figure(figsize=(10,10))
        # plt.ioff() 

        # Loop over all images.
        pbar = ProgressBar()
        for j in pbar(range(len(reduced_files))):
            data = fits.open(reduced_files[j])[0].data
            # Read in some supporting information.
            log_path = pines_path / \
                ('Logs/'+reduced_files[j].name.split('.')[0]+'_log.txt')
            log = pines_log_reader(log_path)
            log_ind = np.where(
                log['Filename'] == reduced_files[j].name.split('_')[0]+'.fits')[0][0]

            header = fits.open(reduced_files[j])[0].header
            date_obs = header['DATE-OBS']
            # Catch a case that can cause datetime strptime to crash; Mimir headers sometimes have DATE-OBS with seconds specified as 010.xx seconds, when it should be 10.xx seconds.
            if len(date_obs.split(':')[-1].split('.')[0]) == 3:
                date_obs = date_obs.split(
                    ':')[0] + ':' + date_obs.split(':')[1] + ':' + date_obs.split(':')[-1][1:]

            if date_obs.split(':')[-1] == '60.00':
                date_obs = date_obs.split(
                    ':')[0]+':'+str(int(date_obs.split(':')[1])+1)+':00.00'
            # Keep a try/except clause here in case other unknown DATE-OBS formats pop up.
            if '.' not in date_obs:
                date_fmt = '%Y-%m-%dT%H:%M:%S'
            else:
                date_fmt = '%Y-%m-%dT%H:%M:%S.%f'
            date = datetime.datetime.strptime(date_obs, date_fmt)

            # Get the closest date master_dark_stddev image for this exposure time.
            # We'll use this to measure read noise and dark current.
            master_dark_stddev = master_dark_stddev_chooser(pines_path/('Calibrations/Darks/Master Darks Stddev/'), header)

            days = date.day + \
                hmsm_to_days(date.hour, date.minute,
                             date.second, date.microsecond)
            jd = date_to_jd(date.year, date.month, days)
            var_df['Filename'][j] = reduced_files[j].name
            var_df['Time UT'][j] = header['DATE-OBS']
            var_df['Time JD UTC'][j] = jd
            var_df['Time BJD TDB'][j] = jd_utc_to_bjd_tdb(
                jd, header['TELRA'], header['TELDEC'])
            var_df['Night Number'][j] = centroided_sources['Night Number'][j]
            var_df['Block Number'][j] = centroided_sources['Block Number'][j]
            var_df['Airmass'][j] = header['AIRMASS']
            var_df['Seeing'][j] = log['X seeing'][np.where(
                log['Filename'] == reduced_files[j].name.split('_')[0]+'.fits')[0][0]]

            # If the shift quality has been flagged, skip this image.
            if log['Shift quality flag'].iloc[log_ind] == 1:
                continue

            # Get the source positions in this image.
            positions = []
            for k in range(len(source_names)):
                positions.append(
                    (centroided_sources[source_names[k]+' Image X'][j], centroided_sources[source_names[k]+' Image Y'][j]))

            # Create an aperture centered on this position with radius (in pixels) of (seeing*multiplicative_factor[j])/plate_scale.
            apertures = CircularAperture(positions, r=(moving_avg_seeing[j]*fact)/plate_scale)

            # Create an annulus centered on this position.
            annuli = CircularAnnulus(positions, r_in=an_in, r_out=an_out)

            photometry_tbl = iraf_style_photometry(
                apertures, annuli, data*gain, master_dark_stddev*gain, header, var_df['Seeing'][j])
            
            # plt.imshow(data, origin='lower', norm=ImageNormalize(data[330:430, 650:750], interval=ZScaleInterval(), stretch=SquaredStretch()))
            # plt.xlim(650, 750)
            # plt.ylim(330, 430)
            # plt.plot(positions[0][0], positions[0][1], 'rx')
            # apertures[0].plot(color='r', lw=2)
            # annuli[0].plot(color='m', lw=2)
            # plt.tight_layout()
            # plt.savefig(output_filename.parent/('save_images/'+str(j).zfill(3)+'.png'), dpi=200)
            # plt.clf()

            for k in range(len(photometry_tbl)):
                var_df[source_names[k]+' Flux'][j] = photometry_tbl['flux'][k]
                var_df[source_names[k] +
                       ' Flux Error'][j] = photometry_tbl['flux_error'][k]
                var_df[source_names[k] +
                       ' Background'][j] = photometry_tbl['background'][k]
                var_df[source_names[k]+' Interpolation Flag'][j] = int(
                    photometry_tbl['interpolation_flag'][k])

        # Write output to file.
        print('Saving multiplicative factor = {} variable aperture photometry output to {}.'.format(
            fact, output_filename))
        print('')
        with open(output_filename, 'w') as f:
            for j in range(len(var_df)):
                # Write in the header.
                if j == 0:
                    f.write('{:>21s}, {:>22s}, {:>17s}, {:>17s}, {:>17s}, {:>17s}, {:>7s}, {:>7s}, '.format(
                        'Filename', 'Time UT', 'Time JD UTC', 'Time BJD TDB', 'Night Number', 'Block Number', 'Airmass', 'Seeing'))
                    for k in range(len(source_names)):
                        if k != len(source_names) - 1:
                            f.write('{:>22s}, {:>28s}, {:>28s}, {:>34s}, '.format(
                                source_names[k]+' Flux', source_names[k]+' Flux Error', source_names[k]+' Background', source_names[k]+' Interpolation Flag'))
                        else:
                            f.write('{:>22s}, {:>28s}, {:>28s}, {:>34s}\n'.format(
                                source_names[k]+' Flux', source_names[k]+' Flux Error', source_names[k]+' Background', source_names[k]+' Interpolation Flag'))

                # Write in Filename, Time UT, Time JD, Airmass, Seeing values.
                format_string = '{:21s}, {:22s}, {:17.9f}, {:17.9f}, {:17d}, {:17d},{:7.2f}, {:7.1f}, '
                # If the seeing value for this image is 'nan' (a string), convert it to a float.
                # TODO: Not sure why it's being read in as a string, fix that.
                if type(var_df['Seeing'][j]) == str:
                    var_df['Seeing'][j] = float(var_df['Seeing'][j])

                # Do a try/except clause for writeout, in case it breaks in the future.
                try:
                    f.write(format_string.format(var_df['Filename'][j], var_df['Time UT'][j], var_df['Time JD UTC'][j], var_df['Time BJD TDB']
                            [j], var_df['Night Number'][j], var_df['Block Number'][j], var_df['Airmass'][j], var_df['Seeing'][j]))
                except:
                    print(
                        'Writeout failed! Inspect quantities you are trying to write out.')
                    breakpoint()

                # Write in Flux, Flux Error, and Background values for every source.
                for i in range(len(source_names)):
                    if i != len(source_names) - 1:
                        format_string = '{:22.5f}, {:28.5f}, {:28.5f}, {:34d}, '
                    else:
                        format_string = '{:22.5f}, {:28.5f}, {:28.5f}, {:34d}\n'
                    try:
                        f.write(format_string.format(var_df[source_names[i]+' Flux'][j], var_df[source_names[i]+' Flux Error']
                                [j], var_df[source_names[i]+' Background'][j], var_df[source_names[i]+' Interpolation Flag'][j]))
                    except:
                        if i != len(source_names) - 1:
                            format_string = '{:22.5f}, {:28.5f}, {:28.5f}, {:34f}, '
                        else:
                            format_string = '{:22.5f}, {:28.5f}, {:28.5f}, {:34f}\n'
                        f.write(format_string.format(var_df[source_names[i]+' Flux'][j], var_df[source_names[i]+' Flux Error']
                                [j], var_df[source_names[i]+' Background'][j], var_df[source_names[i]+' Interpolation Flag'][j]))
        raw_flux_plot(output_filename, mode='night')
        print('')
    return


def basic_psf_phot(short_name, centroided_sources, plots=False):
    warnings.simplefilter("ignore", category=AstropyUserWarning)

    def hmsm_to_days(hour=0, min=0, sec=0, micro=0):
        """
        Convert hours, minutes, seconds, and microseconds to fractional days.

        """
        days = sec + (micro / 1.e6)
        days = min + (days / 60.)
        days = hour + (days / 60.)
        return days / 24.

    def date_to_jd(year, month, day):
        """
        Convert a date to Julian Day.

        Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet', 
            4th ed., Duffet-Smith and Zwart, 2011.

        """
        if month == 1 or month == 2:
            yearp = year - 1
            monthp = month + 12
        else:
            yearp = year
            monthp = month

        # this checks where we are in relation to October 15, 1582, the beginning
        # of the Gregorian calendar.
        if ((year < 1582) or
            (year == 1582 and month < 10) or
                (year == 1582 and month == 10 and day < 15)):
            # before start of Gregorian calendar
            B = 0
        else:
            # after start of Gregorian calendar
            A = math.trunc(yearp / 100.)
            B = 2 - A + math.trunc(A / 4.)

        if yearp < 0:
            C = math.trunc((365.25 * yearp) - 0.75)
        else:
            C = math.trunc(365.25 * yearp)

        D = math.trunc(30.6001 * (monthp + 1))

        jd = B + C + D + day + 1720994.5

        return jd

    def gaussian(p, x, y):
        height, center_x, center_y, width_x, width_y, rotation = p
        rotation = np.deg2rad(rotation)
        x_0 = center_x * np.cos(rotation) - center_y * np.sin(rotation)
        y_0 = center_x * np.sin(rotation) + center_y * np.cos(rotation)

        def rotgauss(x, y):
            xp = x * np.cos(rotation) - y * np.sin(rotation)
            yp = x * np.sin(rotation) + y * np.cos(rotation)
            g = height*np.exp(
                -(((x_0-xp)/width_x)**2 +
                  ((y_0-yp)/width_y)**2)/2.)
            return g

        g = rotgauss(x, y)

        return g

    def moments(data):
        total = np.nansum(data)
        X, Y = np.indices(data.shape)
        center_x = int(np.shape(data)[1]/2)
        center_y = int(np.shape(data)[0]/2)
        row = data[int(center_x), :]
        col = data[:, int(center_y)]
        width_x = np.nansum(np.sqrt(abs((np.arange(col.size)-center_y)**2*col))
                            / np.nansum(col))
        width_y = np.nansum(np.sqrt(abs((np.arange(row.size)-center_x)**2*row))
                            / np.nansum(row))
        height = np.nanmax(data)
        rotation = 0.0
        return height, center_x, center_y, width_x, width_y, rotation

    def errorfunction(p, x, y, data):
        return gaussian(p, x, y) - data

    def fitgaussian(data):
        params = moments(data)
        X, Y = np.indices(data.shape)
        mask = ~np.isnan(data)
        x = X[mask]
        y = Y[mask]
        data = data[mask]
        p, success = optimize.leastsq(errorfunction, params, args=(x, y, data))
        return p

    pines_path = pines_dir_check()
    reduced_path = pines_path/('Objects/'+short_name+'/reduced/')
    reduced_files = np.array(
        natsorted([x for x in reduced_path.glob('*.fits')]))

    centroided_sources.columns = centroided_sources.columns.str.strip()
    source_names = get_source_names(centroided_sources)

    # Declare a new dataframe to hold the information for all targets for this .
    columns = ['Filename', 'Time UT', 'Time JD', 'Airmass', 'Seeing']
    for i in range(0, len(source_names)):
        columns.append(source_names[i]+' Flux')
        columns.append(source_names[i]+' Flux Error')
    psf_df = pd.DataFrame(index=range(len(reduced_files)), columns=columns)
    output_filename = pines_path / \
        ('Objects/'+short_name+'/psf_phot/'+short_name+'_psf_phot.csv')

    for i in range(len(reduced_files)):
        # Read in image data/header.
        file = reduced_files[i]
        data = fits.open(file)[0].data
        header = fits.open(file)[0].header
        print('{}, image {} of {}.'.format(file.name, i+1, len(reduced_files)))

        # Read in some supporting information.
        log_path = pines_path/('Logs/'+file.name.split('.')[0]+'_log.txt')
        log = pines_log_reader(log_path)
        date_obs = header['DATE-OBS']
        # Catch a case that can cause datetime strptime to crash; Mimir headers sometimes have DATE-OBS with seconds specified as 010.xx seconds, when it should be 10.xx seconds.
        if len(date_obs.split(':')[-1].split('.')[0]) == 3:
            date_obs = date_obs.split(
                ':')[0] + ':' + date_obs.split(':')[1] + ':' + date_obs.split(':')[-1][1:]
        # Keep a try/except clause here in case other unknown DATE-OBS formats pop up.
        try:
            date = datetime.datetime.strptime(date_obs, '%Y-%m-%dT%H:%M:%S.%f')
        except:
            print('Header DATE-OBS format does not match the format code in strptime! Inspect/correct the DATE-OBS value.')
            breakpoint()

        days = date.day + \
            hmsm_to_days(date.hour, date.minute, date.second, date.microsecond)
        jd = date_to_jd(date.year, date.month, days)
        psf_df['Filename'][i] = file.name
        psf_df['Time UT'][i] = header['DATE-OBS']
        psf_df['Time JD'][i] = jd
        psf_df['Airmass'][i] = header['AIRMASS']
        psf_df['Seeing'][i] = log['X seeing'][np.where(
            log['Filename'] == file.name.split('_')[0]+'.fits')[0][0]]

        # Read in source centroids for this image
        x = np.zeros(len(source_names))
        y = np.zeros(len(source_names))
        seeing = psf_df['Seeing'][i]

        for j in range(len(source_names)):
            source = source_names[j]
            x[j] = centroided_sources[source+' Image X'][i]
            y[j] = centroided_sources[source+' Image Y'][i]

         # The extract_stars() function requires the input data as an NDData object.
        nddata = NDData(data=data)

        # Create table of good star positions
        stars_tbl = Table()
        stars_tbl['x'] = x
        stars_tbl['y'] = y

        size = 25
        x, y = np.meshgrid(np.arange(0, size), np.arange(0, size))

        # Extract star cutouts.
        stars = extract_stars(nddata, stars_tbl, size=size)

        fitter = fitting.LevMarLSQFitter()

        fig, ax = plt.subplots(nrows=len(stars), ncols=3,
                               sharex=True, sharey=True, figsize=(12, 40))

        # Fit a 2D Gaussian to each star.
        for j in range(len(stars)):
            star = stars[j]
            source = source_names[j]
            mmm_bkg = MMMBackground()
            cutout = star.data - mmm_bkg(star.data)

            # Get the star's centroid position in the cutout.
            dtype = [('x_0', 'f8'), ('y_0', 'f8')]
            pos = Table(data=np.zeros(1, dtype=dtype))
            source_x = stars_tbl['x'][j]
            source_y = stars_tbl['y'][j]
            pos['x_0'] = source_x - int(source_x - size/2 + 1)
            pos['y_0'] = source_y - int(source_y - size/2 + 1)

            parameters = fitgaussian(cutout)
            g2d_fit = gaussian(parameters, x, y)

            avg, med, std = sigma_clipped_stats(cutout)
            im = ax[j, 0].imshow(cutout, origin='lower',
                                 vmin=med-std, vmax=med+8*std)
            divider = make_axes_locatable(ax[j, 0])
            cax = divider.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(im, cax=cax, orientation='vertical')
            ax[j, 0].plot(pos['x_0'], pos['y_0'], 'rx')
            ax[j, 0].set_ylabel(source)
            ax[j, 0].text(pos['x_0'], pos['y_0']+1, '('+str(np.round(source_x, 1)) +
                          ', '+str(np.round(source_y, 1))+')', color='r', ha='center')
            ax[j, 0].axis('off')

            axins = ax[j, 0].inset_axes([0.75, 0.75, 0.25, 0.25])
            axins.set_yticklabels([])
            axins.set_yticks([])
            axins.set_xticklabels([])
            axins.set_xticks([])
            axins.imshow(data, origin='lower', vmin=med-std, vmax=med+8*std)
            axins.plot(source_x, source_y, 'rx')

            im = ax[j, 1].imshow(g2d_fit, origin='lower',
                                 vmin=med-std, vmax=med+8*std)
            divider = make_axes_locatable(ax[j, 1])
            cax = divider.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(im, cax=cax, orientation='vertical')
            ax[j, 1].axis('off')

            avg, med, std = sigma_clipped_stats(cutout - g2d_fit)
            im = ax[j, 2].imshow(
                cutout - g2d_fit, origin='lower', vmin=med-std, vmax=med+8*std)
            divider = make_axes_locatable(ax[j, 2])
            cax = divider.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(im, cax=cax, orientation='vertical')
            ax[j, 2].axis('off')

            if j == 0:
                ax[j, 0].set_title('Data')
                ax[j, 1].set_title('2D Gaussian Model')
                ax[j, 2].set_title('Data - Model')

            plt.tight_layout()
        breakpoint()
        output_filename = pines_path/('Objects/'+short_name+'/basic_psf_phot/' +
                                      reduced_files[i].name.split('_')[0]+'_'+'source_modeling.pdf')
        plt.savefig(output_filename)
        plt.close()

    return


def centroider(short_name, sources, output_plots=False, restore=False, box_w=16, force_output_path='', bin_mins=0.0, shift_tolerance=2.0, time_threshold=0.10):
    """Measures pixel positions of sources in a set of reduced images.

    :param short_name: the target's short name
    :type short_name: str
    :param sources: dataframe of source names and pixel positions, output from ref_star_chooser
    :type sources: pandas DataFrame
    :param output_plots: whether or not to save cutouts of the measured pixel positions, defaults to False
    :type output_plots: bool, optional
    :param restore: whether or not to restore centroid output from a previous run, defaults to False
    :type restore: bool, optional
    :param box_w: the size in pixels of the cutouts used for centroiding, defaults to 16
    :type box_w: int, optional
    :param force_output_path: the top-level path if you don't want to use the default ~/Documents/PINES_analysis_toolkit/ directory for analysis, defaults to ''
    :type force_output_path: str, optional
    :param bin_mins: number of minutes to bin data for staring observations, defaults to 0.0
    :type bin_mins: float, optional
    :param shift_tolerance: number of pixels in x/y a measured source centroid can be from the expected position before alternate methods are tried
    :type shift_tolerance: float, optional
    :return: Saves csv of centroid positions for every source
    :rtype: csv
    """
    # Turn off warnings from Astropy because they're incredibly annoying.
    warnings.simplefilter("ignore", category=AstropyUserWarning)

    #Gaussian filter cutouts for better centroiding SNR
    kernel = Gaussian2DKernel(x_stddev=0.2) 

    matplotlib.use('TkAgg')
    plt.ioff()
    t1 = time.time()

    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()

    kernel = Gaussian2DKernel(x_stddev=1)  # For fixing nans in cutouts.

    # If restore == True, read in existing output and return.
    if restore:
        centroid_df = pd.read_csv(pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv'),
                                  converters={'X Centroids': eval, 'Y Centroids': eval})
        print('Restoring centroider output from {}.'.format(
            pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv')))
        print('')
        return centroid_df

    # Create subdirectories in sources folder to contain output plots.
    if output_plots:
        subdirs = glob(str(pines_path/('Objects/'+short_name+'/sources'))+'/*/')
        # Delete any source directories that are already there.
        for name in subdirs:
            shutil.rmtree(name)

        # Create new source directories.
        for name in sources['Name']:
            source_path = (
                pines_path/('Objects/'+short_name+'/sources/'+name+'/'))
            os.mkdir(source_path)

    # Read in extra shifts, in case the master image wasn't used for source detection.
    extra_shift_path = pines_path / \
        ('Objects/'+short_name+'/sources/extra_shifts.txt')
    extra_shifts = pd.read_csv(extra_shift_path, delimiter=' ', names=[
                               'Extra X shift', 'Extra Y shift'])
    extra_x_shift = extra_shifts['Extra X shift'][0]
    extra_y_shift = extra_shifts['Extra Y shift'][0]

    # Suppress some warnings we don't care about in median combining.
    np.seterr(divide='ignore', invalid='ignore')

    # Get list of reduced files for target.
    reduced_path = pines_path/('Objects/'+short_name+'/reduced')
    reduced_filenames = natsorted(
        [x.name for x in reduced_path.glob('*red.fits')])
    reduced_files = np.array([reduced_path/i for i in reduced_filenames])

    # Declare a new dataframe to hold the centroid information for all sources we want to track.
    columns = []
    columns.append('Filename')
    columns.append('Time JD UTC')
    columns.append('Time BJD TDB')
    columns.append('Night Number')
    columns.append('Block Number')
    columns.append('Seeing')
    columns.append('Airmass')

    # Add x/y positions and cenroid flags for every tracked source
    for i in range(0, len(sources)):
        columns.append(sources['Name'][i]+' Image X')
        columns.append(sources['Name'][i]+' Image Y')
        columns.append(sources['Name'][i]+' Cutout X')
        columns.append(sources['Name'][i]+' Cutout Y')
        columns.append(sources['Name'][i]+' Centroid Warning')

    centroid_df = pd.DataFrame(index=range(
        len(reduced_files)), columns=columns)

    log_path = pines_path/('Logs/')
    log_dates = np.array(
        natsorted([x.name.split('_')[0] for x in log_path.glob('*.txt')]))

    # Make sure we have logs for all the nights of these data. Need them to account for image shifts.
    nights = list(set([i.name.split('.')[0] for i in reduced_files]))
    for i in nights:
        if i not in log_dates:
            print('ERROR: {} not in {}. Download it from the PINES server.'.format(
                i+'_log.txt', log_path))
            breakpoint()

    # Get times before the main loop.
    for j in range(len(reduced_files)):
        file = reduced_files[j]
        header = fits.open(file)[0].header
        time_str = fits.open(file)[0].header['DATE-OBS']

        # Correct some formatting issues that can occur in Mimir time stamps.
        if time_str.split(':')[-1] == '60.00':
            time_str = time_str[0:14] + \
                str(int(time_str.split(':')[-2])+1)+':00.00'
        elif time_str.split(':')[-1] == '010.00':
            time_str = time_str[0:17]+time_str.split(':')[-1][1:]

        if '.' not in time_str:
            time_fmt = '%Y-%m-%dT%H:%M:%S'
        else:
            time_fmt = '%Y-%m-%dT%H:%M:%S.%f'

        jd = julian.to_jd(datetime.datetime.strptime(time_str, time_fmt))
        centroid_df['Time JD UTC'][j] = jd
        centroid_df['Time BJD TDB'][j] = jd_utc_to_bjd_tdb(
            jd, header['TELRA'], header['TELDEC'])

    # From times, generate night and block numbers.
    night_inds = night_splitter(centroid_df['Time BJD TDB'])
    all_night_inds = np.concatenate([np.zeros(
        len(night_inds[i]), dtype='int')+i+1 for i in range(len(night_inds))]).ravel()
    all_block_inds = []
    # On each night...
    for k in range(len(night_inds)):
        # Generate block indices...
        block_inds = block_splitter(centroid_df['Time BJD TDB'][night_inds[k]], bin_mins=bin_mins, time_threshold=time_threshold)
        # Convert them to block numbers...
        night_block_numbers = [
            np.zeros(len(block_inds[i]), dtype='int')+i+1 for i in range(len(block_inds))]
        # And append to the all_block_inds list
        all_block_inds.extend(np.concatenate(night_block_numbers).ravel())
    for j in range(len(reduced_files)):
        centroid_df['Night Number'][j] = all_night_inds[j]
        centroid_df['Block Number'][j] = all_block_inds[j]

    # Measure centroids for each tracked source i in each image j.
    for i in range(len(sources)):
        # Get the initial source position.
        x_pos = sources['Source Detect X'][i]
        y_pos = sources['Source Detect Y'][i]
        print('')
        print('Getting centroids for {}, ({:3.1f}, {:3.1f}) in source detection image. Source {} of {}.'.format(
            sources['Name'][i], x_pos, y_pos, i+1, len(sources)))
        if output_plots:
            print('Saving centroid plots to {}.'.format(pines_path /
                  ('Objects/'+short_name+'/sources/'+sources['Name'][i]+'/')))
        pbar = ProgressBar()
        for j in pbar(range(len(reduced_files))):
            centroid_df[sources['Name'][i]+' Centroid Warning'][j] = 0
            file = reduced_files[j]
            image = fits.open(file)[0].data
            header = fits.open(file)[0].header

            # Get the measured image shift for this image.
            log = pines_log_reader(
                log_path/(file.name.split('.')[0]+'_log.txt'))
            log_ind = np.where(log['Filename'] ==
                            file.name.split('_')[0]+'.fits')[0][0]

            x_shift = float(log['X shift'][log_ind])
            y_shift = float(log['Y shift'][log_ind])

            # Save the filename for readability. Save the seeing for use in variable aperture photometry. Save the time for diagnostic plots.
            if i == 0:
                centroid_df['Filename'][j] = file.name.split('_')[0]+'.fits'
                centroid_df['Seeing'][j] = log['X seeing'][log_ind]
                centroid_df['Airmass'][j] = log['Airmass'][log_ind]

            # Flag indicating if you should not trust the log's shifts. Set to true if x_shift/y_shift are 'nan' or > 30 pixels.
            nan_flag = False

            # If bad shifts were measured for this image, skip.
            if log['Shift quality flag'][log_ind] == 1:
                continue

            if np.isnan(x_shift) or np.isnan(y_shift):
                x_shift = 0
                y_shift = 0
                nan_flag = True

            # If there are clouds, shifts could have been erroneously high...just zero them?
            if abs(x_shift) > 200:
                #x_shift = 0
                nan_flag = True
            if abs(y_shift) > 200:
                #y_shift = 0
                nan_flag = True

            # Apply the shift. NOTE: This relies on having accurate x_shift and y_shift values from the log.
            # If they're incorrect, the cutout will not be in the right place.
            #x_pos = sources['Source Detect X'][i] - x_shift + extra_x_shift
            #y_pos = sources['Source Detect Y'][i] + y_shift - extra_y_shift

            x_pos = sources['Source Detect X'][i] - (x_shift - extra_x_shift)
            y_pos = sources['Source Detect Y'][i] + (y_shift - extra_y_shift)

            # TODO: Make all this its own function.

            # Cutout around the expected position and interpolate over any NaNs (which screw up source detection).
            cutout = interpolate_replace_nans(image[int(y_pos-box_w/2):int(y_pos+box_w/2)+1,int(x_pos-box_w/2):int(x_pos+box_w/2)+1], kernel=Gaussian2DKernel(x_stddev=0.5))
            
            # interpolate_replace_nans struggles with edge pixels, so shave off edge_shave pixels in each direction of the cutout.
            edge_shave = 1
            cutout = cutout[edge_shave:len(cutout)-edge_shave, edge_shave:len(cutout)-edge_shave]

            # #Check for cosmic rays in the cutout. 
            # avg, med, std = sigma_clipped_stats(cutout, sigma=4)
            # counts, bins = np.histogram(cutout.flatten(), bins=20)
            # mode = bins[np.where(counts == max(counts))[0][0]]
            # gap_threshold = 3*std
            # zero_count_bins = np.where(counts == 0)[0]
            # def consecutive(data, stepsize=1):
            #     return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)
            # consecutive_zero_count_bins = consecutive(zero_count_bins)

            # if np.max([len(i) for i in consecutive_zero_count_bins]) > 1:
            #     for v in range(len(consecutive_zero_count_bins)):
            #         start = consecutive_zero_count_bins[v][0]
            #         end = consecutive_zero_count_bins[v][-1]+1
            #         gap = bins[end]-bins[start]
            #         if gap > gap_threshold:
            #             breakpoint()

            #Smooth the cutout with a 2D Gaussian filter for better centroiding SNR.
            cutout = convolve(cutout, kernel, boundary='wrap')

            # Get sigma clipped     stats on the cutout
            vals, lower, upper = sigmaclip(cutout, low=1.5, high=2.5)
            med = np.nanmedian(vals)

            # cs = centroid_sources(cutout-med, (box_w-edge_shave)/2, (box_w-edge_shave)/2)

            # centroid_x_cutout = cs[0][0]
            # centroid_y_cutout = cs[1][0]

            g = centroid_2dg(cutout-med)
            centroid_x_cutout = g[0]
            centroid_y_cutout = g[1]
            
            relative_x_cutout = centroid_x_cutout / (box_w-edge_shave)
            relative_y_cutout = centroid_y_cutout / (box_w-edge_shave)
            
            if (abs(relative_x_cutout - 0.5) > 1) or (abs(relative_y_cutout - 0.5) > 1):
                #breakpoint()
                centroid_x_cutout = (box_w-edge_shave)/2
                centroid_y_cutout = (box_w-edge_shave)/2
                
            # Translate the detected centroid from the cutout coordinates back to the full-frame coordinates.
            centroid_x = centroid_x_cutout + int(x_pos) - box_w/2 + edge_shave
            centroid_y = centroid_y_cutout + int(y_pos) - box_w/2 + edge_shave

            #Make sure you found a centroid in the cutout.
            if (centroid_x_cutout < 0) or (centroid_x_cutout > cutout.shape[1]) or (centroid_y_cutout < 0) or (centroid_y_cutout > cutout.shape[0]):
                print('Centroid found outside of cutout, returning expected position.')
                centroid_x_cutout = cutout.shape[1]/2
                centroid_y_cutout = cutout.shape[0]/2
                centroid_x = x_pos
                centroid_y = y_pos
            
            
            # Check that your measured position is actually on the detector.
            if (centroid_x < 0) or (centroid_y < 0) or (centroid_x > image.shape[1]) or (centroid_y > image.shape[0]):
                # Try a quick mask/interpolation of the cutout.
                mask = make_source_mask(
                    cutout, nsigma=3, npixels=5, dilate_size=3)
                vals, lo, hi = sigmaclip(cutout[~mask])
                bad_locs = np.where((mask == False) & (
                    (cutout > hi) | (cutout < lo)))
                cutout[bad_locs] = np.nan
                cutout = interpolate_replace_nans(
                    cutout, kernel=Gaussian2DKernel(x_stddev=0.5))
                centroid_x, centroid_y = centroid_2dg(cutout - med)
                centroid_x += int(x_pos) - box_w
                centroid_y += int(y_pos) - box_w
                if (centroid_x < 0) or (centroid_y < 0) or (centroid_x > 1023) or (centroid_y > 1023):
                    print(
                        'WARNING: large centroid deviation measured, returning predicted position')
                    print('')
                    centroid_df[sources['Name'][i]+' Centroid Warning'][j] = 1
                    centroid_x = x_pos
                    centroid_y = y_pos
                    # breakpoint()

            # Check to make sure you didn't measure nan's.
            if np.isnan(centroid_x):
                centroid_x = x_pos
                print(
                    'NaN returned from centroid algorithm, defaulting to target position in source_detect_image.')
            if np.isnan(centroid_y):
                centroid_y = y_pos
                print(
                    'NaN returned from centroid algorithm, defaulting to target position in source_detct_image.')

            # Record the image and relative cutout positions.
            centroid_df[sources['Name'][i]+' Image X'][j] = centroid_x
            centroid_df[sources['Name'][i]+' Image Y'][j] = centroid_y
            centroid_df[sources['Name'][i]+' Cutout X'][j] = relative_x_cutout #save coutout as fraction of box_w
            centroid_df[sources['Name'][i]+' Cutout Y'][j] = relative_y_cutout

            if output_plots:
                # Plot
                lock_x = int(centroid_df[sources['Name'][i]+' Image X'][0])
                lock_y = int(centroid_df[sources['Name'][i]+' Image Y'][0])
                norm = ImageNormalize(data=cutout, interval=ZScaleInterval())
                plt.imshow(cutout, origin='lower', norm=norm)
                plt.plot(centroid_x_cutout, centroid_y_cutout, 'rx')
                ap = CircularAperture((centroid_x_cutout, centroid_y_cutout), r=5)
                ap.plot(lw=2, color='b')
                #plt.ylim(lock_y-box_w/2, lock_y+box_w/2-1)
                #plt.xlim(lock_x-box_w/2, lock_x+box_w/2-1)
                plt.title('CENTROID DIAGNOSTIC PLOT\n'+sources['Name'][i]+', '+reduced_files[j].name+' (image '+str(
                    j+1)+' of '+str(len(reduced_files))+')', fontsize=10)
                plt.text(centroid_x_cutout, centroid_y_cutout+0.5, '('+str(np.round(centroid_x_cutout, 1)) +
                         ', '+str(np.round(centroid_y_cutout, 1))+')', color='r', ha='center')
                plot_output_path = (pines_path/('Objects/'+short_name+'/sources/'+sources['Name'][i]+'/'+file.name.split('_')[
                                    0]+'_night_'+str(centroid_df['Night Number'][j]).zfill(2)+'_block_'+str(centroid_df['Block Number'][j]).zfill(2)+'.jpg'))
                plt.gca().set_axis_off()
                plt.subplots_adjust(top=1, bottom=0, right=1,
                                    left=0, hspace=0, wspace=0)
                plt.margins(0, 0)
                plt.gca().xaxis.set_major_locator(plt.NullLocator())
                plt.gca().yaxis.set_major_locator(plt.NullLocator())
                plt.savefig(plot_output_path, bbox_inches='tight',
                            pad_inches=0, dpi=150)
                plt.close()

    output_filename = pines_path / \
        ('Objects/'+short_name+'/sources/target_and_references_centroids.csv')
    # centroid_df.to_csv(pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv'))

    print('Saving centroiding output to {}.'.format(output_filename))
    with open(output_filename, 'w') as f:
        for j in range(len(centroid_df)):
            # Write the header line.
            if j == 0:
                f.write('{:<17s}, '.format('Filename'))
                f.write('{:<15s}, '.format('Time JD UTC'))
                f.write('{:<15s}, '.format('Time BJD TDB'))
                f.write('{:<15s}, '.format('Night Number'))
                f.write('{:<15s}, '.format('Block Number'))
                f.write('{:<6s}, '.format('Seeing'))
                f.write('{:<7s}, '.format('Airmass'))
                for i in range(len(sources['Name'])):
                    n = sources['Name'][i]
                    if i != len(sources['Name']) - 1:
                        f.write('{:<23s}, {:<23s}, {:<24s}, {:<24s}, {:<34s}, '.format(
                            n+' Image X', n+' Image Y', n+' Cutout X', n+' Cutout Y',  n+' Centroid Warning'))
                    else:
                        f.write('{:<23s}, {:<23s}, {:<24s}, {:<24s}, {:<34s}\n'.format(
                            n+' Image X', n+' Image Y', n+' Cutout X', n+' Cutout Y',  n+' Centroid Warning'))

            # Write in the data lines.
            try:
                f.write('{:<17s}, '.format(centroid_df['Filename'][j]))
                f.write('{:<15.7f}, '.format(centroid_df['Time JD UTC'][j]))
                f.write('{:<15.7f}, '.format(centroid_df['Time BJD TDB'][j]))
                f.write('{:<15d}, '.format(centroid_df['Night Number'][j]))
                f.write('{:<15d}, '.format(centroid_df['Block Number'][j]))
                f.write('{:<6.1f}, '.format(float(centroid_df['Seeing'][j])))
                f.write('{:<7.2f}, '.format(centroid_df['Airmass'][j]))
            except:
                breakpoint()

            for i in range(len(sources['Name'])):
                n = sources['Name'][i]
                if i != len(sources['Name']) - 1:
                    format_string = '{:<23.4f}, {:<23.4f}, {:<24.4f}, {:<24.4f}, {:<34d}, '
                else:
                    format_string = '{:<23.4f}, {:<23.4f}, {:<24.4f}, {:<24.4f}, {:<34d}\n'

                f.write(format_string.format(centroid_df[n+' Image X'][j], centroid_df[n+' Image Y'][j],
                        centroid_df[n+' Cutout X'][j], centroid_df[n+' Cutout Y'][j], centroid_df[n+' Centroid Warning'][j]))
    np.seterr(divide='warn', invalid='warn')
    print('')
    print('centroider runtime: {:.2f} minutes.'.format((time.time()-t1)/60))
    print('')
    return centroid_df

def detect_sources(image_path, seeing_fwhm, edge_tolerance, thresh=6.0, plot=False):
    """Finds sources in a Mimir image.

    :param image_path: path to the image
    :type image_path: pathlib.PosixPath
    :param seeing_fwhm: seeing FWHM in arcsec
    :type seeing_fwhm: float
    :param edge_tolerance: how close a source can be to the edge (pixels)
    :type edge_tolerance: float
    :param thresh: significance threshold, defaults to 6.0
    :type thresh: float, optional
    :param plot: whether or not to plot detected sources, defaults to False
    :type plot: bool, optional
    :return: dataframe of sources
    :rtype: pandas DataFrame
    """

    fwhm = seeing_fwhm/0.579  # FIXED

    # Radius of aperture in pixels for doing quick photometry on detected sources.
    ap_rad = 7

    # Read in the image.
    image = fits.open(image_path)[0].data
    header = fits.open(image_path)[0].header

    # Interpolate nans if any in image.
    kernel = Gaussian2DKernel(x_stddev=0.5)
    image = interpolate_replace_nans(image, kernel)

    # Do a quick background model subtraction (makes source detection easier).
    image = bg_2d(image, box_size=32)

    # Get the sigma_clipped_stats for the image.
    avg, med, std = sigma_clipped_stats(image)

    norm = ImageNormalize(
        data=image, interval=ZScaleInterval(), stretch=LinearStretch())

    if plot:
        title = image_path.name
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 9))
        ax.set_aspect('equal')
        im = ax.imshow(image, origin='lower', norm=norm)
        cax = fig.add_axes([0.94, 0.15, 0.015, 0.7])
        fig.colorbar(im, cax=cax, orientation='vertical', label='ADU')
        ax.set_title(title+' Initial Source Detection')

    print('')
    print('Finding sources in {}.'.format(image_path.name))

    # Detect sources using DAOStarFinder.
    daofind = DAOStarFinder(fwhm=fwhm, threshold=thresh*std, sharplo=0.2)

    initial_sources = daofind(image - med)

    # Do a cut based on source sharpness to get rid of some false detections.
    initial_sources.sort('sharpness')
    bad_sharpness_locs = np.where(initial_sources['sharpness'] < 0.3)[0]
    initial_sources.remove_rows(bad_sharpness_locs)

    # Cut sources that are found within edge_tolerance pix of the edges.
    bad_x = np.where((initial_sources['xcentroid'] < edge_tolerance) | (
        initial_sources['xcentroid'] > 1023-edge_tolerance))[0]
    initial_sources.remove_rows(bad_x)
    bad_y = np.where((initial_sources['ycentroid'] < edge_tolerance) | (
        initial_sources['ycentroid'] > 1023-edge_tolerance))[0]
    initial_sources.remove_rows(bad_y)

    # Cut sources near y = 512, these are frequently bad.
    bad_512 = np.where((initial_sources['ycentroid'] > 506) & (
        initial_sources['ycentroid'] < 518))
    initial_sources.remove_rows(bad_512)

    # Cut sources near the top and bottom edges to avoid issues with the Mimir "ski jump" feature that can sometimes occur.
    bad_ski_jump = np.where((initial_sources['ycentroid'] < 100) | (
        initial_sources['ycentroid'] > 924))
    initial_sources.remove_rows(bad_ski_jump)

    # Do quick photometry on the remaining sources.
    positions = [(initial_sources['xcentroid'][i], initial_sources['ycentroid'][i])
                 for i in range(len(initial_sources))]
    apertures = CircularAperture(positions, r=ap_rad)
    phot_table = aperture_photometry(image-med, apertures)

    # Cut based on brightness.
    phot_table.sort('aperture_sum')
    cutoff = 1*std*np.pi*ap_rad**2
    bad_source_locs = np.where(phot_table['aperture_sum'] < cutoff)
    phot_table.remove_rows(bad_source_locs)
    initial_sources.remove_rows(bad_source_locs)

    if plot:
        # Plot detected sources.
        # TODO: indicate saturated sources, sources near edge, etc. with different color markers.
        ax.plot(phot_table['xcenter'], phot_table['ycenter'],
                'ro', markerfacecolor='none')

    #print('Found {} sources.'.format(len(phot_table)))
    # Resort remaining sources so that the brightest are listed firsts.
    sources = phot_table[::-1].to_pandas()
    return sources

def epsf_phot(target, centroided_sources, plots=False):
    def hmsm_to_days(hour=0, min=0, sec=0, micro=0):
        """
        Convert hours, minutes, seconds, and microseconds to fractional days.

        """
        days = sec + (micro / 1.e6)
        days = min + (days / 60.)
        days = hour + (days / 60.)
        return days / 24.

    def date_to_jd(year, month, day):
        """
        Convert a date to Julian Day.

        Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet', 
            4th ed., Duffet-Smith and Zwart, 2011.

        """
        if month == 1 or month == 2:
            yearp = year - 1
            monthp = month + 12
        else:
            yearp = year
            monthp = month

        # this checks where we are in relation to October 15, 1582, the beginning
        # of the Gregorian calendar.
        if ((year < 1582) or
            (year == 1582 and month < 10) or
                (year == 1582 and month == 10 and day < 15)):
            # before start of Gregorian calendar
            B = 0
        else:
            # after start of Gregorian calendar
            A = math.trunc(yearp / 100.)
            B = 2 - A + math.trunc(A / 4.)

        if yearp < 0:
            C = math.trunc((365.25 * yearp) - 0.75)
        else:
            C = math.trunc(365.25 * yearp)

        D = math.trunc(30.6001 * (monthp + 1))

        jd = B + C + D + day + 1720994.5

        return jd

    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    reduced_path = pines_path/('Objects/'+short_name+'/reduced/')
    reduced_filenames = natsort.natsorted(
        [x.name for x in reduced_path.glob('*.fits')])
    reduced_files = np.array([reduced_path/i for i in reduced_filenames])

    centroided_sources.columns = centroided_sources.columns.str.strip()
    source_names = natsort.natsorted(list(set([i.split(' ')[0]+' '+i.split(
        ' ')[1] for i in centroided_sources.keys() if (i[0] == '2') or (i[0] == 'R')])))

    # Create output plot directories for each source.
    if plots:
        for name in source_names:
            # If the folders are already there, delete them.
            source_path = (
                pines_path/('Objects/'+short_name+'/psf_phot/'+name+'/'))
            if source_path.exists():
                shutil.rmtree(source_path)
            # Create folders.
            os.mkdir(source_path)

    # Declare a new dataframe to hold the information for all targets for this .
    columns = ['Filename', 'Time UT', 'Time JD', 'Airmass', 'Seeing']
    for i in range(0, len(source_names)):
        columns.append(source_names[i]+' Flux')
        columns.append(source_names[i]+' Flux Error')
    psf_df = pd.DataFrame(index=range(len(reduced_files)), columns=columns)
    output_filename = pines_path / \
        ('Objects/'+short_name+'/psf_phot/'+short_name+'_psf_phot.csv')

    for i in range(0, len(reduced_files)):
        # Read in image data/header.
        file = reduced_files[i]
        data = fits.open(file)[0].data
        header = fits.open(file)[0].header
        print('{}, image {} of {}.'.format(file.name, i+1, len(reduced_files)))

        # Read in some supporting information.
        log_path = pines_path/('Logs/'+file.name.split('.')[0]+'_log.txt')
        log = pines_log_reader(log_path)
        date_obs = header['DATE-OBS']
        # Catch a case that can cause datetime strptime to crash; Mimir headers sometimes have DATE-OBS with seconds specified as 010.xx seconds, when it should be 10.xx seconds.
        if len(date_obs.split(':')[-1].split('.')[0]) == 3:
            date_obs = date_obs.split(
                ':')[0] + ':' + date_obs.split(':')[1] + ':' + date_obs.split(':')[-1][1:]
        date = datetime.datetime.strptime(date_obs, '%Y-%m-%dT%H:%M:%S.%f')

        days = date.day + \
            hmsm_to_days(date.hour, date.minute, date.second, date.microsecond)
        jd = date_to_jd(date.year, date.month, days)
        psf_df['Filename'][i] = file.name
        psf_df['Time UT'][i] = header['DATE-OBS']
        psf_df['Time JD'][i] = jd
        psf_df['Airmass'][i] = header['AIRMASS']
        psf_df['Seeing'][i] = log['X seeing'][np.where(
            log['Filename'] == file.name.split('_')[0]+'.fits')[0][0]]

        # Read in source centroids for this image
        x = np.zeros(len(source_names))
        y = np.zeros(len(source_names))
        for j in range(len(source_names)):
            source = source_names[j]
            x[j] = centroided_sources[source+' X'][i]
            y[j] = centroided_sources[source+' Y'][i]

        # Extract pixel cutouts of our stars, so let’s explicitly exclude stars that are too close to the image boundaries (because they cannot be extracted).
        size = 13
        hsize = (size - 1) / 2
        #mask = ((x > hsize) & (x < (data.shape[1] -1 - hsize)) & (y > hsize) & (y < (data.shape[0] -1 - hsize)) & (y > 100) & (y < 923))

        # Create table of good star positions
        stars_tbl = Table()
        stars_tbl['x'] = x
        stars_tbl['y'] = y

        # Subtract background (star cutouts from which we build the ePSF must have background subtracted).
        mean_val, median_val, std_val = sigma_clipped_stats(data, sigma=2.)
        data -= median_val

        # Replace nans in data using Gaussian.
        # kernel = Gaussian2DKernel(x_stddev=0.5)
        # data = interpolate_replace_nans(data, kernel)

        # The extract_stars() function requires the input data as an NDData object.
        nddata = NDData(data=data)

        # Extract star cutouts.
        stars = extract_stars(nddata, stars_tbl, size=size)

        # Plot.
        nrows = 5
        ncols = 5
        fig, ax = plt.subplots(nrows=nrows, ncols=ncols,
                               figsize=(10, 10), squeeze=True)
        ax = ax.ravel()
        for j in range(len(stars)):
            norm = simple_norm(stars[j], 'log', percent=99.)
            ax[j].imshow(stars[j].data, norm=norm,
                         origin='lower', cmap='viridis')

        pdb.set_trace()

        # Construct the ePSF using the star cutouts.
        epsf_fitter = EPSFFitter()
        epsf_builder = EPSFBuilder(
            maxiters=4, progress_bar=False, fitter=epsf_fitter)

        try:
            epsf, fitted_stars = epsf_builder(stars)
            output_filename = pines_path / \
                ('Objects/'+short_name+'/psf_phot/'+short_name+'_psf_phot.csv')

            for j in range(len(stars)):
                star = stars[j]
                source_name = source_names[j]
                sigma_psf = 1.85

                dtype = [('x_0', 'f8'), ('y_0', 'f8')]
                pos = Table(data=np.zeros(1, dtype=dtype))
                source_x = stars_tbl['x'][j]
                source_y = stars_tbl['y'][j]
                pos['x_0'] = source_x - int(source_x - size/2 + 1)
                pos['y_0'] = source_y - int(source_y - size/2 + 1)

                daogroup = DAOGroup(4.0*sigma_psf*gaussian_sigma_to_fwhm)
                mmm_bkg = MMMBackground()
                photometry = BasicPSFPhotometry(group_maker=daogroup,
                                                bkg_estimator=mmm_bkg,
                                                psf_model=epsf,
                                                fitter=LevMarLSQFitter(),
                                                fitshape=(13, 13),
                                                aperture_radius=4.)

                result_tab = photometry(image=star, init_guesses=pos)
                residual_image = photometry.get_residual_image()
                psf_df[source_name+' Flux'][i] = result_tab['flux_fit'][0]
                psf_df[source_name+' Flux Error'][i] = result_tab['flux_unc'][0]

                if plots:
                    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(12, 4))
                    im = ax[0].imshow(star, origin='lower')
                    divider = make_axes_locatable(ax[0])
                    cax = divider.append_axes('right', size='5%', pad=0.05)
                    fig.colorbar(im, cax=cax, orientation='vertical')
                    ax[0].plot(result_tab['x_fit'][0],
                               result_tab['y_fit'][0], 'rx')
                    ax[0].set_title('Data')

                    im2 = ax[1].imshow(epsf.data, origin='lower')
                    ax[1].set_title('EPSF Model')
                    divider = make_axes_locatable(ax[1])
                    cax = divider.append_axes('right', size='5%', pad=0.05)
                    fig.colorbar(im2, cax=cax, orientation='vertical')

                    im3 = ax[2].imshow(residual_image, origin='lower')
                    ax[2].set_title('Residual Image')
                    divider = make_axes_locatable(ax[2])
                    cax = divider.append_axes('right', size='5%', pad=0.05)
                    fig.colorbar(im3, cax=cax, orientation='vertical')
                    plt.suptitle(
                        source_name+'\n'+reduced_files[i].name+', image '+str(i+1)+' of '+str(len(reduced_files)))
                    plt.subplots_adjust(wspace=0.5, top=0.95, bottom=0.05)
                    plot_output_name = pines_path / \
                        ('Objects/'+short_name+'/psf_phot/' +
                         source_name+'/'+str(i).zfill(4)+'.jpg')
                    plt.savefig(plot_output_name)
                    plt.close()
        except:
            print('')
            print('EPSF BUILDER FAILED, SKIPPING IMAGE.')
            print('')
        # Plot the ePSF.
        # plt.figure()
        # norm = simple_norm(epsf.data, 'log', percent=99.)
        # plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
        # cb = plt.colorbar()
        # plt.tight_layout()

    print('Saving psf photometry output to {}.'.format(output_filename))
    with open(output_filename, 'w') as f:
        for j in range(len(psf_df)):
            if j == 0:
                f.write('{:>21s}, {:>22s}, {:>17s}, {:>7s}, {:>7s}, '.format(
                    'Filename', 'Time UT', 'Time JD', 'Airmass', 'Seeing'))
                for i in range(len(source_names)):
                    if i != len(source_names) - 1:
                        f.write('{:>20s}, {:>26s}, '.format(
                            source_names[i]+' Flux', source_names[i]+' Flux Error'))
                    else:
                        f.write('{:>20s}, {:>26s}\n'.format(
                            source_names[i]+' Flux', source_names[i]+' Flux Error'))

            format_string = '{:21s}, {:22s}, {:17.9f}, {:7.2f}, {:7.1f}, '

            # If the seeing value for this image is 'nan' (a string), convert it to a float.
            # TODO: Not sure why it's being read in as a string, fix that.
            if type(psf_df['Seeing'][j]) == str:
                psf_df['Seeing'][j] = float(psf_df['Seeing'][j])

            # Do a try/except clause for writeout, in case it breaks in the future.
            try:
                f.write(format_string.format(psf_df['Filename'][j], psf_df['Time UT']
                        [j], psf_df['Time JD'][j], psf_df['Airmass'][j], psf_df['Seeing'][j]))
            except:
                print('Writeout failed! Inspect quantities you are trying to write out.')
                pdb.set_trace()
            for i in range(len(source_names)):
                if i != len(source_names) - 1:
                    format_string = '{:20.11f}, {:26.11f}, '
                else:
                    format_string = '{:20.11f}, {:26.11f}\n'

                f.write(format_string.format(
                    psf_df[source_names[i]+' Flux'][j], psf_df[source_names[i]+' Flux Error'][j]))
    print('')
    return

def ref_star_chooser(short_name, source_detect_image, guess_position=(700., 382.), radius_check=6., non_linear_limit=3300.,
                     dimness_tolerance=0.5, brightness_tolerance=10., closeness_tolerance=12., distance_from_target=900., edge_tolerance=50., exclude_lower_left=False, restore=False,
                     source_detect_plot=False, force_output_path=''):
    """Chooses suitable reference stars for a target in a specified source_detect_image.

    :param short_name: the short name for the target
    :type short_name: str
    :param source_detect_image: name of the reduced image in which you want to find reference stars
    :type source_detect_image: str
    :param guess_position: guess position for the target in the source_detect_image, defaults to (700.,382.)
    :type guess_position: tuple, optional
    :param radius_check: the radius in pixels used to perform photometry to compare target and reference brightnesses, defaults to 6.
    :type radius_check: float, optional
    :param non_linear_limit: ADU value above which references are considered to be in the non-linear limit of the detecto, defaults to 3300.
    :type non_linear_limit: float, optional
    :param dimness_tolerance: minimum multiple of the target's measured brightness that a reference is allowed to have, defaults to 0.5
    :type dimness_tolerance: float, optional
    :param brightness_tolerance: maximum multiple of the target's measured brightness that a reference is allowed to have, defaults to 10.
    :type brightness_tolerance: float, optional
    :param closeness_tolerance: the closest distance in pixels that a reference star can be to another detected source and still be considered as a reference, defaults to 12.
    :type closeness_tolerance: float, optional
    :param distance_from_target: the furthest distance in pixels that a reference star can be from the target and still be considered, defaults to 900.
    :type distance_from_target: float, optional
    :param edge_tolerance: the closest distance in pixels that a reference can be to the edge of the detector and still be considered, defaults to 50.
    :type edge_tolerance: float, optional
    :param exclude_lower_left: whether or not to exclude reference stars from the lower left quadrant (due to occasional Mimir 'bars' issue), defaults to False
    :type exclude_lower_left: bool, optional
    :param restore: whether or not to restore references from previous output, defaults to False
    :type restore: bool, optional
    :param source_detect_plot: whenther or not to plot all detected sources, defaults to False
    :type source_detect_plot: bool, optional
    :param force_output_path: top-level 'force output path', used if you want to use a folder other than ~/Documents/PINES_analysis_toolkit/, defaults to ''
    :type force_output_path: str, optional
    :raises ValueError: If the measured seeing in your selected source_detect_image is Nan
    :raises ValueError: If the measured x/y shift in your selected_source_detect_image is > 20 pixels
    :return: Saves a plot of the target/detected reference stars to the object's 'sources' directory. Saves .csv of target/reference pixel positions in the object's 'sources' directory.
    :rtype: plot/csv
    """

    plt.ion()
    # Get your local PINES directory
    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()

    if restore:
        output_df = pd.read_csv(
            pines_path/('Objects/'+short_name+'/sources/target_and_references_source_detection.csv'))
        print('')
        print('Restoring ref_star_chooser output from {}.'.format(
            pines_path/('Objects/'+short_name+'/sources/target_and_references_source_detection.csv')))
        print('')
        return output_df

    # Find reduced files in the directory.
    data_dir = pines_path/('Objects/'+short_name+'/reduced/')
    reduced_files = np.array(natsorted([x for x in data_dir.glob('*.fits')]))

    # If the directory doesn't exist, or there are no reduced files, return nothing.
    if (not data_dir.exists()) or (len(reduced_files) == 0):
        print('ERROR: No reduced images exist for {}.'.format(short_name))
        return

    # Set path to source directory.
    source_dir = pines_path/('Objects/'+short_name+'/sources/')

    source_detect_image_ind = np.where(
        [i.name == source_detect_image for i in reduced_files])[0][0]
    source_detect_image_path = reduced_files[source_detect_image_ind]
    source_frame = source_detect_image_path.name.split('_')[0]+'.fits'
    log_name = source_frame.split('.')[0]+'_log.txt'
    log = pines_log_reader(pines_path/('Logs/'+log_name))
    source_frame_ind = np.where(log['Filename'] == source_frame)[0][0]
    extra_x_shift = log['X shift'][source_frame_ind]
    extra_y_shift = log['Y shift'][source_frame_ind]
    source_detect_seeing = float(log['X seeing'][source_frame_ind])

    # Make sure the seeing value isn't a NaN.
    if np.isnan(source_detect_seeing):
        raise ValueError(
            'Seeing value = NaN in log. Try a different source_detect_image_ind in call to ref_star_chooser.')

    if (float(extra_x_shift) > 20) or (float(extra_y_shift) > 20):
        raise ValueError(
            'Measured x or y shift > 20 pixels. Try a different source_detect_image_ind in call to ref_star_chooser.')

    extra_shift_path = (source_dir/'extra_shifts.txt')
    with open(extra_shift_path, 'w') as file:
        file.write(str(extra_x_shift)+' '+str(extra_y_shift))

    #source_detect_seeing = 3.5

    # Detect sources in the image.
    sources = detect_sources(source_detect_image_path, source_detect_seeing,
                             edge_tolerance, plot=source_detect_plot, thresh=3.0)

    # Identify the target in the image using guess_position.
    target_id = target_finder(sources, guess_position)

    # Now determine ids of suitable reference stars.
    image = fits.open(source_detect_image_path)[0].data

    # Interpolate nans if any in image.
    kernel = Gaussian2DKernel(x_stddev=0.5)
    image = interpolate_replace_nans(image, kernel)

    suitable_ref_ids = []
    bad_ref_ids = [target_id]
    # Get a value for the brightness of the target with a radius_check aperture. We don't want reference stars dimmer than dimness_tolerance * targ_flux_estimates.
    target_ap = CircularAperture(
        (sources['xcenter'][target_id], sources['ycenter'][target_id]), r=radius_check)
    target_an = CircularAnnulus(
        (sources['xcenter'][target_id], sources['ycenter'][target_id]), r_in=12, r_out=30)
    mask = target_an.to_mask(method='center')
    an_data = mask.multiply(image)
    an_data_1d = an_data[an_data != 0]
    vals, lower, upper = sigmaclip(an_data_1d)
    bg_estimate = np.median(vals)
    targ_flux_estimate = aperture_photometry(
        image, target_ap)['aperture_sum'] - bg_estimate * target_ap.area

    for i in range(len(sources)):
        if i != target_id:
            potential_ref_loc = (sources['xcenter'][i], sources['ycenter'][i])
            ap = CircularAperture(potential_ref_loc, r=radius_check)
            an = CircularAnnulus(potential_ref_loc, r_in=12, r_out=30)
            mask = an.to_mask(method='center')
            an_data = mask.multiply(image)
            an_data_1d = an_data[an_data != 0]
            vals, lower, upper = sigmaclip(an_data_1d)
            bg_estimate = np.median(vals)

            pixels = ap.to_mask()
            sub_image = pixels.cutout(image)

            dists = []
            for j in range(len(sources)):
                dists.append(np.sqrt((sources['xcenter'][i]-sources['xcenter'][j])**2+(
                    sources['ycenter'][i]-sources['ycenter'][j])**2))
            dists = np.array(dists)

            # Check if potential ref is near non-linear regime.
            non_linear_flag = len(
                np.where(sub_image > non_linear_limit)[0]) == 0
            # Check if potential ref is bright enough.
            dimness_flag = ((aperture_photometry(image, ap)[
                            'aperture_sum'] - bg_estimate * ap.area) > dimness_tolerance*targ_flux_estimate)[0]
            brightness_flag = ((aperture_photometry(image, ap)[
                               'aperture_sum'] - bg_estimate * ap.area) < brightness_tolerance*targ_flux_estimate)[0]
            closeness_flag = (
                len(np.where(dists[np.where(dists != 0)] < closeness_tolerance)[0]) == 0)
            proximity_flag = (np.sqrt((sources['xcenter'][target_id]-sources['xcenter'][i])**2+(
                sources['ycenter'][target_id]-sources['ycenter'][i])**2) < distance_from_target)
            if exclude_lower_left:
                if (sources['xcenter'][i] < 512) & (sources['ycenter'][i] < 512):
                    lower_left_flag = False
                else:
                    lower_left_flag = True
            else:
                lower_left_flag = True

            if (non_linear_flag) & (dimness_flag) & (brightness_flag) & (closeness_flag) & (proximity_flag) & (lower_left_flag):
                suitable_ref_ids.append(i)
            else:
                bad_ref_ids.append(i)

    #ax.plot(sources['xcenter'][suitable_ref_ids], sources['ycenter'][suitable_ref_ids],'yx')
    print('Found {} suitable reference stars.'.format(len(suitable_ref_ids)))
    print('')

    if len(suitable_ref_ids) < 2:
        print('Not enough reference stars found. Try loosening reference star criteria.')
        return

    target = sources.iloc[[target_id]].reset_index(
        drop=True).drop(columns=['id'])
    suitable_refs = sources.drop(bad_ref_ids).reset_index(
        drop=True).drop(columns=['id'])

    output_df = pd.DataFrame(
        columns=['Name', 'Source Detect X', 'Source Detect Y'])
    output_df = output_df.append(
        {'Name': short_name, 'Source Detect X': sources['xcenter'][target_id], 'Source Detect Y': sources['ycenter'][target_id]}, ignore_index=True)
    for i in range(len(suitable_refs)):
        output_df = output_df.append({'Name': 'Reference '+str(
            i+1), 'Source Detect X': sources['xcenter'][suitable_ref_ids[i]], 'Source Detect Y': sources['ycenter'][suitable_ref_ids[i]]}, ignore_index=True)

    stats = sigma_clipped_stats(image)
    fig, ax = plt.subplots(figsize=(10, 9))
    ax.set_aspect('equal')
    im = ax.imshow(
        image, vmin=stats[1], vmax=stats[1]+7*stats[2], origin='lower', cmap='Greys')
    cax = fig.add_axes([0.87, 0.15, 0.035, 0.7])
    fig.colorbar(im, cax=cax, orientation='vertical', label='ADU')
    ax.set_title(source_detect_image_path.name)
    ax.plot(output_df['Source Detect X'][0], output_df['Source Detect Y'][0],
            marker='o', color='m', mew=2, ms=12, ls='', mfc='None', label='Target')
    # Plot selected references
    for i in range(1, len(output_df)):
        ax.plot(output_df['Source Detect X'][i], output_df['Source Detect Y'][i],
                marker='o', color='b', mew=2, ms=12, ls='', mfc='None', label='Reference')
        ax.text(output_df['Source Detect X'][i]+7,
                output_df['Source Detect Y'][i]+5, str(i), fontsize=14, color='r')
    # Plot detected sources
    for i in range(len(sources)):
        ax.plot(sources['xcenter'], sources['ycenter'], 'rx',
                label='Detected sources', alpha=0.6, ms=3)

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.1, 1.1))
    plt.show()

    # Sometimes, the source detection returns things that are clearly not reference stars.
    # Allow the user to remove them here.
    ans = input(
        'Enter IDs of references to remove, separated by commas (e.g.: 1,4,8,14).\nIf none to remove, hit enter:  ')
    if ans != '':
        ans = ans.split(',')
        ref_ids_to_remove = []
        for i in range(len(ans)):
            ref_ids_to_remove.append(int(ans[i]))

        # Reverse sort the list.
        ref_ids_to_remove = np.sort(ref_ids_to_remove)[::-1]
        # Drop those references from output_df
        for i in ref_ids_to_remove:
            output_df = output_df.drop(index=i)
        output_df = output_df.reset_index(drop=True)
        # Rename the remaining references.
        for i in range(1, len(output_df)):
            name = output_df['Name'][i]
            if int(name.split(' ')[1]) != i:
                name = 'Reference '+str(i)
                output_df['Name'][i] = name
        # Replot everything.
        ax.cla()
        im = ax.imshow(image, origin='lower',
                       vmin=stats[1], vmax=stats[1]+7*stats[2], cmap='Greys')
        cax = fig.add_axes([0.87, 0.15, 0.035, 0.7])
        fig.colorbar(im, cax=cax, orientation='vertical', label='ADU')
        ax.set_title(source_detect_image_path.name)
        ax.plot(output_df['Source Detect X'][0], output_df['Source Detect Y'][0],
                marker='o', color='m', mew=2, ms=12, ls='', mfc='None', label='Target')
        for i in range(1, len(output_df)):
            ax.plot(output_df['Source Detect X'][i], output_df['Source Detect Y'][i],
                    marker='o', color='b', mew=2, ms=12, ls='', mfc='None', label='Reference')
            ax.text(output_df['Source Detect X'][i]+7,
                    output_df['Source Detect Y'][i]+5, str(i), color='r', fontsize=14)
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(),
                  bbox_to_anchor=(1.1, 1.1))
    plt.show()

    ans = input('Happy with reference star selection? y/n: ')
    if ans == 'y':
        print('')
        print('Saving target/reference info to {}.'.format(source_dir /
              'target_and_references_source_detection.csv'))
        output_df.to_csv(
            source_dir/'target_and_references_source_detection.csv', index=False)
        print('')
        print('Saving target and references image to {}.'.format(
            source_dir/('target_and_refs.png')))
        plt.savefig(source_dir/('target_and_refs.png'))
        plt.close('all')
        return output_df
    elif ans == 'n':
        raise ValueError(
            'Try changing arguments of ref_star_chooser or detect_sources and try again.')
    else:
        return

def target_finder(sources, guess_position):
    """Identifies the ID of the target given an initial guess of its pixel position.

    :param sources: dataframe of sources 
    :type sources: pandas DataFrame
    :param guess_position: tuple of (x,y) pixel positions for the guess location
    :type guess_position: tuple of floats
    :return: target ID
    :rtype: int
    """
    distances = np.sqrt((guess_position[0] - sources['xcenter'])
                        ** 2 + (guess_position[1] - sources['ycenter'])**2)
    target_id = np.where(distances == np.min(distances))[0][0]
    return target_id

def twoD_Gaussian(x, y, sigma_x, sigma_y, rho, x_cen, y_cen, amplitude, offset):
        
    covar = np.zeros((2,2))
    covar[0,0] = sigma_x**2
    covar[1,1] = sigma_y**2
    covar[0,1] = rho*sigma_x*sigma_y
    covar[1,0] = covar[0,1]
    inv_covar = np.linalg.inv(covar)

    ln_Gaussian = (x - x_cen)**2 * inv_covar[0,0] + 2*(x - x_cen)*(y - y_cen)*inv_covar[0,1] + (y - y_cen)**2 * inv_covar[1,1]
    
    return offset + amplitude * np.exp(-ln_Gaussian/2.0)

def joint_psf_fitter(cutouts, gain=8.21, plate_scale=0.579):
    '''Does a 2D gaussian fit jointly to a set of image cutouts (n_stars x nxpix x nypix)
    '''
    #Inner fit fits a single star with sigma and rho fixed
    def inner_min_function(inner_params,x,y,cutout,err,sigma_x,sigma_y,rho):
        x_cen=inner_params[0]
        y_cen=inner_params[1]
        amplitude=inner_params[2]
        offset=inner_params[3]
        model = twoD_Gaussian(x, y, sigma_x, sigma_y, rho, x_cen, y_cen, amplitude, offset)
        good = np.isfinite(cutout.flatten())
        return (cutout.flatten()[good]-model.flatten()[good])/err.flatten()[good]

    #Outer fit does sigma_x, sigma_y and rho, calling inner_fit for each star
    def outer_min_function(outer_params,xs,ys,nstars,cutouts,errs):
        sigma_x = outer_params[0]
        sigma_y = outer_params[1]
        rho = outer_params[2]
        inner_ini_guess = np.zeros(4)
        for i in np.arange(nstars):
            inner_ini_guess[0] = np.mean(xs[i,:,:])
            inner_ini_guess[1] = np.mean(ys[i,:,:])
            inner_ini_guess[2] = np.nanmax(cutouts[i,:,:]) - np.nanmin(cutouts[i,:,:])
            inner_ini_guess[3] = np.nanmin(cutouts[i,:,:])
            inner_res_lsq = least_squares(inner_min_function,
                                        inner_ini_guess,
                                        args=(xs[i,:,:],
                                        ys[i,:,:],
                                        cutouts[i,:,:],
                                        errs[i,:,:],
                                        sigma_x,
                                        sigma_y,
                                        rho))
            if i==0:
                residuals = np.array(inner_res_lsq['fun']).flatten()
            else:
                residuals = np.concatenate((residuals,np.array(inner_res_lsq['fun']).flatten()))

        return residuals

    def joint_PSF_fit(cutouts, errs, xs, ys, plate_scale):
        
        nstars = cutouts.shape[0]

        outer_ini_guess = np.zeros(3)
        lower_bounds = np.zeros(3)
        upper_bounds = np.zeros(3)

        #x and y sigma guesses and bounds
        outer_ini_guess[0:2] = 3/plate_scale/2.355 #sigma corresponds to 2.5" seeing
        lower_bounds[0:2] = 0.5/plate_scale/2.355 #0.5 arcseconds bounds on seeing
        upper_bounds[0:2] = 10/plate_scale/2.355 #10 arcseconds boudsn on seeing

        #rho guess and bounds
        outer_ini_guess[2] = 0.1#assume no rotation to the Gaussian yet (rho=0)
        lower_bounds[2] = -0.99999999
        upper_bounds[2] = 0.99999999

        #make the bounds a "2 tuple" so least_squares is happy
        bounds = (lower_bounds,upper_bounds)

        #do the fit for sigmas and rho
        res_lsq = least_squares(outer_min_function, outer_ini_guess, args=(xs, ys, nstars, cutouts, errs), diff_step=0.001, verbose=0)

        #now, re-run fit for each star to get the individual params
        sigma_x=res_lsq.x[0]
        sigma_y=res_lsq.x[1]
        rho = res_lsq.x[2]
        xcens = np.zeros(nstars)
        ycens = np.zeros(nstars)
        offsets = np.zeros(nstars)
        amplitudes = np.zeros(nstars)

        inner_ini_guess = np.zeros(4)

        for i in np.arange(nstars):
            inner_ini_guess[0] = np.mean(xs[i,:,:])
            inner_ini_guess[1] = np.mean(ys[i,:,:])
            inner_ini_guess[2] = np.nanmax(cutouts[i,:,:]) - np.nanmin(cutouts[i,:,:])
            inner_ini_guess[3] = np.nanmin(cutouts[i,:,:])
            res_lsq = least_squares(inner_min_function,
                                    inner_ini_guess,
                                    args=(xs[i,:,:],
                                        ys[i,:,:],
                                        cutouts[i,:,:],
                                        errs[i,:,:],
                                        sigma_x,
                                        sigma_y,
                                        rho),
                                    verbose=0)
            xcens[i] = res_lsq.x[0]
            ycens[i] = res_lsq.x[1]
            amplitudes[i] = res_lsq.x[2]
            offsets[i] = res_lsq.x[3]

        return sigma_x,sigma_y,rho,xcens,ycens,amplitudes,offsets

    #Make an error matrix to match cutouts. Assume ph noise:
    errs = np.sqrt(cutouts*gain)/gain


    #Get sizes
    nxpix = cutouts.shape[1]
    nypix = cutouts.shape[2]

    #No xs and ys provided, so just them all on the the same grid
    x = np.resize(np.arange(nxpix),(nxpix,nypix))
    y = np.transpose(np.resize(np.arange(nypix),(nypix,nxpix)))

    #resize to make a 3d array that matches cutouts
    xs = np.resize(x,cutouts.shape)
    ys = np.resize(y,cutouts.shape)

    #run the code
    sigma_x,sigma_y,rho,xcens,ycens,amplitudes,offsets = joint_PSF_fit(cutouts, errs, xs, ys, plate_scale)

    fit_params = {'sigma_x':sigma_x, 'sigma_y':sigma_y, 'rho':rho, 'xcens':xcens, 'ycens':ycens, 'amplitudes':amplitudes, 'offsets':offsets}

    return fit_params, xs, ys

def cutout_array_creator(image, x_c, y_c, cutout_w=11):
    cutout_array = np.zeros((len(x_c), 2*cutout_w, 2*cutout_w))
    for i in range(len(x_c)):
        cutout = image[int(y_c[i]-cutout_w):int(y_c[i]+cutout_w), int(x_c[i]-cutout_w):int(x_c[i]+cutout_w)]
        cutout_array[i] = cutout
    return cutout_array
