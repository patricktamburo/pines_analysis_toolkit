import pdb
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader
from pines_analysis_toolkit.utils.quick_plot import quick_plot
import numpy as np
import natsort
from photutils import CircularAperture, CircularAnnulus, aperture_photometry, make_source_mask
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
import datetime 
import math
from photutils.utils import calc_total_error
from astropy.table import Table, QTable
import astropy.units as u
from scipy.stats import sigmaclip
import shutil
import os
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from progressbar import ProgressBar
from astropy.visualization import ZScaleInterval, ImageNormalize, SquaredStretch, MinMaxInterval
from pines_analysis_toolkit.data.master_dark_stddev_chooser import master_dark_stddev_chooser
from pines_analysis_toolkit.utils.quick_plot import quick_plot as qp 
from astropy.modeling import models, fitting
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pines_analysis_toolkit.utils.get_source_names import get_source_names
import matplotlib
matplotlib.use('TkAgg')

def hmsm_to_days(hour=0,min=0,sec=0,micro=0):
    """
    Convert hours, minutes, seconds, and microseconds to fractional days.
    
    """
    days = sec + (micro / 1.e6)
    days = min + (days / 60.)
    days = hour + (days / 60.)
    return days / 24.

def date_to_jd(year,month,day):
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

def iraf_style_photometry(phot_apertures, bg_apertures, data, dark_std_data, header, seeing, bg_method='mean', epadu=1.0):
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
        raise ValueError('Invalid background method, choose either mean, median, or mode')
    
    #Create a list to hold the flux for each source.
    aperture_sum = []
    interpolation_flags = np.zeros(len(phot_apertures.positions), dtype='bool')

    for i in range(len(phot_apertures.positions)):
        pos = phot_apertures.positions[i]
        #Cutout around the source position
        cutout_w = 15
        x_pos = pos[0]
        y_pos = pos[1]
        cutout = data[int((y_pos-cutout_w)):int(y_pos+cutout_w)+1, int(x_pos-cutout_w):int(x_pos+cutout_w)+1]
        x_cutout = x_pos - np.floor(x_pos - cutout_w)
        y_cutout = y_pos - np.floor(y_pos - cutout_w) 
        ap = CircularAperture((x_cutout, y_cutout), r=phot_apertures.r)

        #Cut out the pixels JUST inside the aperture, and check if there are NaNs there. If so, interpolate over NaNs in the cutout. 
        ap_mask = ap.to_mask(method='exact')
        ap_cut = ap_mask.cutout(cutout)

        bad_sum = np.sum(np.isnan(ap_cut))

        if bad_sum > 0:
            bads = np.where(np.isnan(cutout))
            bad_dists = np.sqrt((bads[0] - y_cutout)**2 + (bads[1] - x_cutout)**2)

            #Check if any bad pixels fall within the aperture. If so, set the interpolation flag to True for this source. 
            if np.sum(bad_dists < phot_apertures.r+1):

                # if np.sum(bad_dists < 1) == 0:
                #     #ONLY interpolate if bad pixels lay away from centroid position by at least a pixel. 
                #     #2D gaussian fitting approach
                #     #Set up a 2D Gaussian model to interpolate the bad pixel values in the cutout. 
                #     model_init = models.Const2D(amplitude=np.nanmedian(cutout))+models.Gaussian2D(amplitude=np.nanmax(cutout), x_mean=x_cutout, y_mean=y_cutout, x_stddev=seeing, y_stddev=seeing)
                #     xx, yy = np.indices(cutout.shape) #2D grids of x and y coordinates
                #     mask = ~np.isnan(cutout) #Find locations where the cutout has *good* values (i.e. not NaNs). 
                #     x = xx[mask] #Only use coordinates at these good values.
                #     y = yy[mask]
                #     cutout_1d = cutout[mask] #Make a 1D cutout using only the good values. 
                #     fitter = fitting.LevMarLSQFitter()
                #     model_fit = fitter(model_init, x, y, cutout_1d) #Fit the model to the 1d cutout.
                #     cutout[~mask] = model_fit(xx,yy)[~mask] #Interpolate the pixels in the cutout using the 2D Gaussian fit. 
                #     pdb.set_trace()
                #     #TODO: interpolate_replace_nans with 2DGaussianKernel probably gives better estimation of *background* pixels. 
                # else:
                #     interpolation_flags[i] = True

                #2D gaussian fitting approach
                #Set up a 2D Gaussian model to interpolate the bad pixel values in the cutout. 
                model_init = models.Const2D(amplitude=np.nanmedian(cutout))+models.Gaussian2D(amplitude=np.nanmax(cutout), x_mean=x_cutout, y_mean=y_cutout, x_stddev=seeing, y_stddev=seeing)
                xx, yy = np.indices(cutout.shape) #2D grids of x and y coordinates
                mask = ~np.isnan(cutout) #Find locations where the cutout has *good* values (i.e. not NaNs). 
                x = xx[mask] #Only use coordinates at these good values.
                y = yy[mask]
                cutout_1d = cutout[mask] #Make a 1D cutout using only the good values. 
                fitter = fitting.LevMarLSQFitter()
                model_fit = fitter(model_init, x, y, cutout_1d) #Fit the model to the 1d cutout.

                # #Uncomment this block to show inerpolation plots.
                # norm = ImageNormalize(cutout, interval=ZScaleInterval())
                # plt.ion()
                # fig, ax = plt.subplots(1, 4, figsize=(10,4), sharex=True, sharey=True)
                # ax[0].imshow(cutout, origin='lower', norm=norm)
                # ax[0].set_title('Data')
                # ax[1].imshow(model_fit(xx,yy), origin='lower', norm=norm)
                # ax[1].set_title('2D Gaussian Model')

                cutout[~mask] = model_fit(xx,yy)[~mask] #Interpolate the pixels in the cutout using the 2D Gaussian fit. 

                # ax[2].imshow(cutout, origin='lower', norm=norm)
                # ax[2].set_title('Data w/ Bad\nPixels Replaced')
                # ax[3].imshow(cutout-model_fit(xx,yy), origin='lower')
                # ax[3].set_title('Residuals')
                # pdb.set_trace()

                interpolation_flags[i] = True
                
                # #Gaussian convolution approach.
                # cutout = interpolate_replace_nans(cutout, kernel=Gaussian2DKernel(x_stddev=0.5))


        phot_source = aperture_photometry(cutout, ap)
        # if np.isnan(phot_source['aperture_sum'][0]):
        #     pdb.set_trace()

        aperture_sum.append(phot_source['aperture_sum'][0])

    #Add positions/fluxes to a table
    xcenter = phot_apertures.positions[:,0]*u.pix
    ycenter = phot_apertures.positions[:,1]*u.pix
    phot = QTable([xcenter, ycenter, aperture_sum], names=('xcenter', 'ycenter', 'aperture_sum'))

    #Now measure the background around each source. 
    mask = make_source_mask(data, nsigma=3, npixels=5, dilate_size=7) #Make a mask to block out any sources that might show up in the annuli and bias them.
    bg_phot = aperture_stats_tbl(~mask*data, bg_apertures, sigma_clip=True) #Pass the data with sources masked out to the bg calculator. 
    ap_area = phot_apertures.area
    
    bg_method_name = 'aperture_{}'.format(bg_method)
    background = bg_phot[bg_method_name]
    flux = phot['aperture_sum'] - background * ap_area

    # Need to use variance of the sources for Poisson noise term in error computation.
    flux_error = compute_phot_error(flux, bg_phot, bg_method, ap_area, exptime, dark_std_data, phot_apertures, epadu)

    mag = -2.5 * np.log10(flux)
    mag_err = 1.0857 * flux_error / flux

    # Make the final table
    X, Y = phot_apertures.positions.T
    stacked = np.stack([X, Y, flux, flux_error, mag, mag_err, background, interpolation_flags], axis=1)
    names = ['X', 'Y', 'flux', 'flux_error', 'mag', 'mag_error', 'background', 'interpolation_flag']

    final_tbl = Table(data=stacked, names=names)

    #Check for nans
    if sum(np.isnan(final_tbl['flux'])) > 0:
        bad_locs = np.where(np.isnan(final_tbl['flux']))[0]
        #pdb.set_trace()

    return final_tbl
    
def compute_phot_error(flux, bg_phot, bg_method, ap_area, exptime, dark_std_data, phot_apertures, epadu=1.0):
    """Computes the flux errors using the DAOPHOT style computation.
        Includes photon noise from the source, background terms, dark current, and read noise."""

    #See eqn. 1 from Broeg et al. (2005): https://ui.adsabs.harvard.edu/abs/2005AN....326..134B/abstract
    flux_variance_term = flux / epadu #This is just the flux in photons. 
    bg_variance_term_1 = ap_area * (bg_phot['aperture_std'])**2
    bg_variance_term_2 = (ap_area * bg_phot['aperture_std'])**2 / bg_phot['aperture_area']

    #Measure combined read noise + dark current using the dark standard deviation image. 
    #TODO: Check that this is correct. 
    dark_rn_term = np.zeros(len(flux_variance_term))
    for i in range(len(phot_apertures)):
        ap = phot_apertures[i]
        dark_rn_ap = ap.to_mask().multiply(dark_std_data)
        dark_rn_ap = dark_rn_ap[dark_rn_ap != 0]
        dark_rn_term[i] = ap_area * (np.median(dark_rn_ap)**2)
    
    variance = flux_variance_term + bg_variance_term_1 + bg_variance_term_2 + dark_rn_term
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
    names = ['X','Y','aperture_mean','aperture_median','aperture_mode','aperture_std', 'aperture_area']
    
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
        values = cutout[cutout != 0] #Unwrap the annulus into a 1D array. 
        values = values[~np.isnan(values)] #Ignore any NaNs in the annulus. 
        if sigma_clip:
            values, clow, chigh = sigmaclip(values, low=3, high=3) #Sigma clip the 1D annulus values. 
        mean = np.nanmean(values)
        median = np.nanmedian(values)
        std = np.nanstd(values)
        mode = 3 * median - 2 * mean
        actual_area = (~np.isnan(values)).sum()
        return (mean, median, mode, std, actual_area)
    

def fixed_aper_phot(target, centroided_sources, ap_radii, an_in=12., an_out=30., plots=False, gain=8.21, qe=0.9):
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
        plots (bool, optional): Whether or not to output surface plots. Images output to aper_phot directory within the object directory.
        gain (float, optional): The gain of the detector in e-/ADU.
        qe (float, optional): The quantum efficiency of the detector.
    Outputs:
        Saves aperture photometry csv to PINES_analysis_toolkit/Objects/short_name/aper_phot/ for each aperture.
	TODO:
    '''

    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    
    #Remove any leading/trailing spaces in the column names. 
    centroided_sources.columns = centroided_sources.columns.str.lstrip()
    centroided_sources.columns = centroided_sources.columns.str.rstrip()

    #Get list of reduced files for target. 
    reduced_path = pines_path/('Objects/'+short_name+'/reduced')
    reduced_filenames = natsort.natsorted([x.name for x in reduced_path.glob('*red.fits')])
    reduced_files = np.array([reduced_path/i for i in reduced_filenames])
    
    #source_names = natsort.natsorted(list(set([i.replace('X','').replace('Y','').replace('Centroid Warning','').strip() for i in centroided_sources.keys() if i != 'Filename'])))
    source_names = get_source_names(centroided_sources)
    
    #Create output plot directories for each source.
    if plots:
        #Camera angles for surface plots
        azim_angles = np.linspace(0, 360*1.5, len(reduced_files)) % 360
        elev_angles = np.zeros(len(azim_angles)) + 25
        for name in source_names:
            #If the folders are already there, delete them. 
            source_path = (pines_path/('Objects/'+short_name+'/aper_phot/'+name+'/'))
            if source_path.exists():
                shutil.rmtree(source_path)
            #Create folders.
            os.mkdir(source_path)

    #Loop over all aperture radii. 
    for ap in ap_radii:
        print('Doing fixed aperture photometry for {}, aperture radius = {:1.1f} pix, inner annulus radius = {} pix, outer annulus radius = {} pix.'.format(target, ap, an_in, an_out))

        #Declare a new dataframe to hold the information for all targets for this aperture. 
        columns = ['Filename', 'Time UT', 'Time JD', 'Airmass', 'Seeing']
        for i in range(0, len(source_names)):
            columns.append(source_names[i]+' Flux')
            columns.append(source_names[i]+' Flux Error')
            columns.append(source_names[i]+' Background')
            columns.append(source_names[i]+' Interpolation Flag')

        ap_df = pd.DataFrame(index=range(len(reduced_files)), columns=columns)
        output_filename = pines_path/('Objects/'+short_name+'/aper_phot/'+short_name+'_fixed_aper_phot_{:1.1f}_pix_radius.csv'.format(float(ap)))

        #Loop over all images.
        pbar = ProgressBar()
        for j in pbar(range(len(reduced_files))):
            data = fits.open(reduced_files[j])[0].data

            #Read in some supporting information.
            log_path = pines_path/('Logs/'+reduced_files[j].name.split('.')[0]+'_log.txt')
            log = pines_log_reader(log_path)
            log_ind = np.where(log['Filename'] == reduced_files[j].name.split('_')[0]+'.fits')[0][0]

            

            header = fits.open(reduced_files[j])[0].header
            date_obs = header['DATE-OBS']
            #Catch a case that can cause datetime strptime to crash; Mimir headers sometimes have DATE-OBS with seconds specified as 010.xx seconds, when it should be 10.xx seconds. 
            if len(date_obs.split(':')[-1].split('.')[0]) == 3:
                date_obs = date_obs.split(':')[0] + ':' + date_obs.split(':')[1] + ':' + date_obs.split(':')[-1][1:]
            
            if date_obs.split(':')[-1] == '60.00':
                date_obs = date_obs.split(':')[0]+':'+str(int(date_obs.split(':')[1])+1)+':00.00'
            #Keep a try/except clause here in case other unknown DATE-OBS formats pop up. 
            try:
                date = datetime.datetime.strptime(date_obs, '%Y-%m-%dT%H:%M:%S.%f')
            except:
                print('Header DATE-OBS format does not match the format code in strptime! Inspect/correct the DATE-OBS value.')
                pdb.set_trace()
            
            #Get the closest date master_dark_stddev image for this exposure time.
            #We'll use this to measure read noise and dark current. 
            date_str = date_obs.split('T')[0].replace('-','')
            master_dark_stddev = master_dark_stddev_chooser(pines_path/('Calibrations/Darks/Master Darks Stddev/'), header)
            
            days = date.day + hmsm_to_days(date.hour,date.minute,date.second,date.microsecond)
            jd = date_to_jd(date.year,date.month,days)
            ap_df['Filename'][j] = reduced_files[j].name
            ap_df['Time UT'][j] = header['DATE-OBS']
            ap_df['Time JD'][j] = jd
            ap_df['Airmass'][j] = header['AIRMASS']
            ap_df['Seeing'][j] = log['X seeing'][log_ind]
            
            #If the shift quality has been flagged, skip this image. 
            if log['Shift quality flag'].iloc[log_ind] == 1:
                continue

            #Get the source positions in this image.
            positions = []
            for i in range(len(source_names)):
                positions.append((float(centroided_sources[source_names[i]+' Image X'][j]), float(centroided_sources[source_names[i]+' Image Y'][j])))

            #Create an aperture centered on this position with radius = ap. 
            try:
                apertures = CircularAperture(positions, r=ap)
            except:
                pdb.set_trace()

            #Create an annulus centered on this position. 
            annuli = CircularAnnulus(positions, r_in=an_in, r_out=an_out)

            photometry_tbl = iraf_style_photometry(apertures, annuli, data*gain, master_dark_stddev*gain, header, ap_df['Seeing'][j])

            for i in range(len(photometry_tbl)):
                ap_df[source_names[i]+' Flux'][j] = photometry_tbl['flux'][i] 
                ap_df[source_names[i]+' Flux Error'][j] = photometry_tbl['flux_error'][i]
                ap_df[source_names[i]+' Background'][j] = photometry_tbl['background'][i]
                ap_df[source_names[i]+' Interpolation Flag'][j] = int(photometry_tbl['interpolation_flag'][i])

            #Make surface plots.
            if plots:
                for i in range(len(photometry_tbl)):
                    x_p = photometry_tbl['X'][i]
                    y_p = photometry_tbl['Y'][i]

                    fig = plt.figure()
                    ax = fig.add_subplot(111, projection='3d')
                    xx, yy = np.meshgrid(np.arange(int(x_p)-10, int(x_p)+10+1), np.arange(int(y_p)-10, int(y_p)+10+1))
                    theta = np.linspace(0, 2 * np.pi, 201)
                    y_circ = ap*np.cos(theta)+y_p
                    x_circ = ap*np.sin(theta)+x_p
                    vmin = np.nanmedian(data[yy,xx])
                    vmax = vmin + 2.5*np.nanstd(data[yy,xx])
                    ax.plot_surface(xx, yy, data[yy,xx], cmap=cm.viridis, alpha=0.8, rstride=1, cstride=1, edgecolor='k', lw=0.2, vmin=vmin, vmax=vmax)
                    current_z = ax.get_zlim()
                    ax.set_zlim(current_z[0]-150, current_z[1])
                    current_z = ax.get_zlim()
                    cset = ax.contourf(xx, yy, data[yy,xx], zdir='z', offset=current_z[0], cmap=cm.viridis)
                    ax.plot(x_circ, y_circ, np.zeros(len(x_circ))+current_z[0], color='r', lw=2, zorder=100)
                    ax.set_xlabel('X')
                    ax.set_ylabel('Y')
                    ax.set_zlabel('Counts')

                    ax.set_title('SURFACE DIAGNOSTIC PLOT, '+', Ap. = '+str(ap)+'\n'+source_names[i]+', '+reduced_files[j].name+' (image '+str(j+1)+' of '+str(len(reduced_files))+')')
                    ax.view_init(elev=elev_angles[j], azim=azim_angles[j])
                    plot_output_path = (pines_path/('Objects/'+short_name+'/aper_phot/'+source_names[i]+'/'+str(j).zfill(4)+'.jpg'))
                    plt.tight_layout()
                    plt.savefig(plot_output_path)
                    plt.close()
           
        #Write output to file. 
        print('Saving ap = {:1.1f} aperture photometry output to {}.'.format(ap,output_filename))
        print('')
        with open(output_filename, 'w') as f:
            for j in range(len(ap_df)):
                #Write in the header. 
                if j == 0:
                    f.write('{:>21s}, {:>22s}, {:>17s}, {:>7s}, {:>7s}, '.format('Filename', 'Time UT', 'Time JD', 'Airmass', 'Seeing'))
                    for i in range(len(source_names)):
                        if i != len(source_names) - 1:
                            f.write('{:>22s}, {:>28s}, {:>28s}, {:>34s}, '.format(source_names[i]+' Flux', source_names[i]+' Flux Error', source_names[i]+' Background', source_names[i]+' Interpolation Flag'))
                        else:
                            f.write('{:>22s}, {:>28s}, {:>28s}, {:>34s}\n'.format(source_names[i]+' Flux', source_names[i]+' Flux Error', source_names[i]+' Background', source_names[i]+' Interpolation Flag'))


                #Write in Filename, Time UT, Time JD, Airmass, Seeing values.
                format_string = '{:21s}, {:22s}, {:17.9f}, {:7.2f}, {:7.1f}, '
                #If the seeing value for this image is 'nan' (a string), convert it to a float. 
                #TODO: Not sure why it's being read in as a string, fix that. 
                if type(ap_df['Seeing'][j]) == str:
                    ap_df['Seeing'][j] = float(ap_df['Seeing'][j])

                #Do a try/except clause for writeout, in case it breaks in the future. 
                try:
                    f.write(format_string.format(ap_df['Filename'][j], ap_df['Time UT'][j], ap_df['Time JD'][j], ap_df['Airmass'][j], ap_df['Seeing'][j]))
                except:
                    print('Writeout failed! Inspect quantities you are trying to write out.')
                    pdb.set_trace()

                #Write in Flux, Flux Error, and Background values for every source. 
                for i in range(len(source_names)):                    
                    if i != len(source_names) - 1:
                        format_string = '{:22.5f}, {:28.5f}, {:28.5f}, {:34d}, '
                    else:
                        format_string = '{:22.5f}, {:28.5f}, {:28.5f}, {:34d}\n'
                    try:
                        f.write(format_string.format(ap_df[source_names[i]+' Flux'][j], ap_df[source_names[i]+' Flux Error'][j], ap_df[source_names[i]+' Background'][j], ap_df[source_names[i]+' Interpolation Flag'][j]))
                    except:
                        if i != len(source_names) - 1:
                            format_string = '{:22.5f}, {:28.5f}, {:28.5f}, {:34f}, '
                        else:
                            format_string = '{:22.5f}, {:28.5f}, {:28.5f}, {:34f}\n'
                        f.write(format_string.format(ap_df[source_names[i]+' Flux'][j], ap_df[source_names[i]+' Flux Error'][j], ap_df[source_names[i]+' Background'][j], ap_df[source_names[i]+' Interpolation Flag'][j]))

    print('')
    return 

def variable_aper_phot(target, centroided_sources, multiplicative_factors, an_in=12., an_out=30., plots=False, gain=8.21, qe=0.9, plate_scale=0.579):
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    
    #Remove any leading/trailing spaces in the column names. 
    centroided_sources.columns = centroided_sources.columns.str.lstrip()
    centroided_sources.columns = centroided_sources.columns.str.rstrip()

    #Get list of reduced files for target. 
    reduced_path = pines_path/('Objects/'+short_name+'/reduced')
    reduced_filenames = natsort.natsorted([x.name for x in reduced_path.glob('*red.fits')])
    reduced_files = np.array([reduced_path/i for i in reduced_filenames])
    
    #Get source names. 
    source_names = get_source_names(centroided_sources)
    
    #Get seeing.
    seeing = np.array(centroided_sources['Seeing'])

    #Loop over multiplicative factors 
    for i in range(len(multiplicative_factors)):
        fact = multiplicative_factors[i]
        print('Doing variable aperture photometry for {}, multiplicative seeing factor = {}, inner annulus radius = {} pix, outer annulus radius = {} pix.'.format(target, fact, an_in, an_out))

        #Declare a new dataframe to hold the information for all targets for this aperture. 
        columns = ['Filename', 'Time UT', 'Time JD', 'Airmass', 'Seeing']
        for j in range(0, len(source_names)):
            columns.append(source_names[j]+' Flux')
            columns.append(source_names[j]+' Flux Error')
            columns.append(source_names[j]+' Background')
            columns.append(source_names[j]+' Interpolation Flag')

        var_df = pd.DataFrame(index=range(len(reduced_files)), columns=columns)
        output_filename = pines_path/('Objects/'+short_name+'/aper_phot/'+short_name+'_variable_aper_phot_'+str(float(fact))+'_seeing_factor.csv')

        #Loop over all images.
        pbar = ProgressBar()
        for j in pbar(range(len(reduced_files))):
            data = fits.open(reduced_files[j])[0].data
            #Read in some supporting information.
            log_path = pines_path/('Logs/'+reduced_files[j].name.split('.')[0]+'_log.txt')
            log = pines_log_reader(log_path)
            header = fits.open(reduced_files[j])[0].header
            date_obs = header['DATE-OBS']
            #Catch a case that can cause datetime strptime to crash; Mimir headers sometimes have DATE-OBS with seconds specified as 010.xx seconds, when it should be 10.xx seconds. 
            if len(date_obs.split(':')[-1].split('.')[0]) == 3:
                date_obs = date_obs.split(':')[0] + ':' + date_obs.split(':')[1] + ':' + date_obs.split(':')[-1][1:]
            
            if date_obs.split(':')[-1] == '60.00':
                date_obs = date_obs.split(':')[0]+':'+str(int(date_obs.split(':')[1])+1)+':00.00'
            #Keep a try/except clause here in case other unknown DATE-OBS formats pop up. 
            try:
                date = datetime.datetime.strptime(date_obs, '%Y-%m-%dT%H:%M:%S.%f')
            except:
                print('Header DATE-OBS format does not match the format code in strptime! Inspect/correct the DATE-OBS value.')
                pdb.set_trace()
            
            #Get the closest date master_dark_stddev image for this exposure time.
            #We'll use this to measure read noise and dark current. 
            date_str = date_obs.split('T')[0].replace('-','')
            master_dark_stddev = master_dark_stddev_chooser(pines_path/('Calibrations/Darks/Master Darks Stddev/'), header)
        
            days = date.day + hmsm_to_days(date.hour,date.minute,date.second,date.microsecond)
            jd = date_to_jd(date.year,date.month,days)
            var_df['Filename'][j] = reduced_files[j].name
            var_df['Time UT'][j] = header['DATE-OBS']
            var_df['Time JD'][j] = jd
            var_df['Airmass'][j] = header['AIRMASS']
            var_df['Seeing'][j] = log['X seeing'][np.where(log['Filename'] == reduced_files[j].name.split('_')[0]+'.fits')[0][0]]
            
            #Get the source positions in this image.
            positions = []
            for k in range(len(source_names)):
                positions.append((centroided_sources[source_names[k]+' Image X'][j], centroided_sources[source_names[k]+' Image Y'][j]))

            #Create an aperture centered on this position with radius (in pixels) of (seeing*multiplicative_factor[j])/plate_scale. 
            apertures = CircularAperture(positions, r=(seeing[j]*fact)/plate_scale)

            #Create an annulus centered on this position. 
            annuli = CircularAnnulus(positions, r_in=an_in, r_out=an_out)

            photometry_tbl = iraf_style_photometry(apertures, annuli, data*gain, master_dark_stddev*gain, header, var_df['Seeing'][j])

            for k in range(len(photometry_tbl)):
                var_df[source_names[k]+' Flux'][j] = photometry_tbl['flux'][k] 
                var_df[source_names[k]+' Flux Error'][j] = photometry_tbl['flux_error'][k]
                var_df[source_names[k]+' Background'][j] = photometry_tbl['background'][k]
                var_df[source_names[k]+' Interpolation Flag'][j] = int(photometry_tbl['interpolation_flag'][k])
        
        #Write output to file. 
        print('Saving multiplicative factor = {} variable aperture photometry output to {}.'.format(fact,output_filename))
        print('')
        with open(output_filename, 'w') as f:
            for j in range(len(var_df)):
                #Write in the header. 
                if j == 0:
                    f.write('{:>21s}, {:>22s}, {:>17s}, {:>7s}, {:>7s}, '.format('Filename', 'Time UT', 'Time JD', 'Airmass', 'Seeing'))
                    for k in range(len(source_names)):
                        if k != len(source_names) - 1:
                            f.write('{:>22s}, {:>28s}, {:>28s}, {:>34s}, '.format(source_names[k]+' Flux', source_names[k]+' Flux Error', source_names[k]+' Background', source_names[k]+' Interpolation Flag'))
                        else:
                            f.write('{:>22s}, {:>28s}, {:>28s}, {:>34s}\n'.format(source_names[k]+' Flux', source_names[k]+' Flux Error', source_names[k]+' Background', source_names[k]+' Interpolation Flag'))


                #Write in Filename, Time UT, Time JD, Airmass, Seeing values.
                format_string = '{:21s}, {:22s}, {:17.9f}, {:7.2f}, {:7.1f}, '
                #If the seeing value for this image is 'nan' (a string), convert it to a float. 
                #TODO: Not sure why it's being read in as a string, fix that. 
                if type(var_df['Seeing'][j]) == str:
                    var_df['Seeing'][j] = float(var_df['Seeing'][j])

                #Do a try/except clause for writeout, in case it breaks in the future. 
                try:
                    f.write(format_string.format(var_df['Filename'][j], var_df['Time UT'][j], var_df['Time JD'][j], var_df['Airmass'][j], var_df['Seeing'][j]))
                except:
                    print('Writeout failed! Inspect quantities you are trying to write out.')
                    pdb.set_trace()

                #Write in Flux, Flux Error, and Background values for every source. 
                for k in range(len(source_names)):                    
                    if k != len(source_names) - 1:
                        format_string = '{:22.5f}, {:28.5f}, {:28.5f}, {:34d}, '
                    else:
                        format_string = '{:22.5f}, {:28.5f}, {:28.5f}, {:34d}\n'
                    f.write(format_string.format(var_df[source_names[k]+' Flux'][j], var_df[source_names[k]+' Flux Error'][j], var_df[source_names[k]+' Background'][j], var_df[source_names[k]+' Interpolation Flag'][j]))
        print('')
    return