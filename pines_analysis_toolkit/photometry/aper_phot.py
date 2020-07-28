import pdb
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader
from pines_analysis_toolkit.utils.quick_plot import quick_plot
import numpy as np
import natsort
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
import datetime 
import math
from photutils.utils import calc_total_error
from astropy.table import Table
from scipy.stats import sigmaclip
# from .background_median import aperture_stats_tbl

'''Authors:
		Patrick Tamburo, Boston University, June 2020
	Purpose:
        Performs aperture photometry on a set of reduced images given dataframe of source positions.
        The iraf_style_photometry, compute_phot_error, perture_stats_tbl, and calc_aperture_mmm routines are from Varun Bajaj on github:
            https://github.com/spacetelescope/wfc3_photometry/blob/master/photometry_tools/photometry_with_errors.py. 
	Inputs:
        target (str): The target's full 2MASS name.
        sources (pandas dataframe): List of source names, x and y positions in every image. 
        ap_radii (list of floats): List of aperture radii in pixels for which aperture photometry wil be performed. 
        an_in (float, optional): The inner radius of the annulus used to estimate background, in pixels. 
        an_out (float, optional): The outer radius of the annulus used to estimate background, in pixels. 
    Outputs:
        Saves aperture photometry csv to PINES_analysis_toolkit/Objects/short_name/aper_phot/ for each aperture.
	TODO:
        Make sure Varun's routines are working properly. Why is photometric uncertainty < sqrt(flux)? 
'''

def aper_phot(target, centroided_sources, ap_radii, an_in=12., an_out=30.):
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
    
    def iraf_style_photometry(phot_apertures, bg_apertures, data, error_array=None, bg_method='mode', epadu=1.0):
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
        error_array: array, optional
            The array of pixelwise error of the data.  If none, the
            Poisson noise term in the error computation will just be the
            square root of the flux/epadu. If not none, the
            aperture_sum_err column output by aperture_photometry
            (divided by epadu) will be used as the Poisson noise term.
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

        if bg_method not in ['mean', 'median', 'mode']:
            raise ValueError('Invalid background method, choose either mean, median, or mode')
        
        phot = aperture_photometry(data, phot_apertures, error=error_array)
        bg_phot = aperture_stats_tbl(data, bg_apertures, sigma_clip=True)
        ap_area = phot_apertures.area
        bg_method_name = 'aperture_{}'.format(bg_method)

        flux = phot['aperture_sum'] - bg_phot[bg_method_name] * ap_area

        # Need to use variance of the sources for Poisson noise term in error computation.
        #This means error needs to be squared. If no error_array error = flux ** .5
        if error_array is not None:
            flux_error = compute_phot_error(phot['aperture_sum_err']**2.0, bg_phot, bg_method, ap_area, epadu)
        else:
            flux_error = compute_phot_error(flux, bg_phot, bg_method, ap_area, epadu)

        mag = -2.5 * np.log10(flux)
        mag_err = 1.0857 * flux_error / flux

        # Make the final table
        X, Y = phot_apertures.positions.T
        stacked = np.stack([X, Y, flux, flux_error, mag, mag_err], axis=1)
        names = ['X', 'Y', 'flux', 'flux_error', 'mag', 'mag_error']

        final_tbl = Table(data=stacked, names=names)

        #Check for nans
        if sum(np.isnan(final_tbl['flux'])) > 0:
            bad_locs = np.where(np.isnan(final_tbl['flux']))[0]
            quick_plot(data)
            plt.plot(final_tbl['X'][bad_locs], final_tbl['Y'][bad_locs], 'rx')

            pdb.set_trace()

        return final_tbl
        
    def compute_phot_error(flux_variance, bg_phot, bg_method, ap_area, epadu=1.0):
        """Computes the flux errors using the DAOPHOT style computation"""
        #See eqn. 1 from Broeg et al. (2005): https://ui.adsabs.harvard.edu/abs/2005AN....326..134B/abstract
        flux_variance_term = flux_variance / epadu
        bg_variance_term_1 = ap_area * (bg_phot['aperture_std'])**2
        bg_variance_term_2 = (ap_area * bg_phot['aperture_std'])**2 / bg_phot['aperture_area']
        variance = flux_variance_term + bg_variance_term_1 + bg_variance_term_2
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
            values = cutout[cutout != 0]
            if sigma_clip:
                values, clow, chigh = sigmaclip(values, low=3, high=3)
            mean = np.mean(values)
            median = np.median(values)
            std = np.std(values)

            mode = 3 * median - 2 * mean
            actual_area = (~np.isnan(values)).sum()
            return (mean, median, mode, std, actual_area)
    
    plt.ion()
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    
    #Remove any leading/trailing spaces in the column names. 
    centroided_sources.columns = centroided_sources.columns.str.lstrip()
    centroided_sources.columns = centroided_sources.columns.str.rstrip()

    #Get list of reduced files for target. 
    reduced_path = pines_path/('Objects/'+short_name+'/reduced')
    reduced_files = np.array(natsort.natsorted([x for x in reduced_path.glob('*.fits')]))

    source_names = natsort.natsorted(list(set([i[0:-2].replace('X','').replace('Y','').rstrip().lstrip() for i in centroided_sources.keys()])))
    #Old way of doing the above line. 
    # for i in centroided_sources.keys():
    #     name = i.split(' ')[0]+' '+i.split(' ')[1]
    #     pdb.set_trace()
    #     if name not in source_names:
    #         source_names.append(name)
    
    if len(source_names) != len(centroided_sources.keys())/2:
        print('ERROR: something going wrong grabbing source names. ')
        return
        
    #Loop over all aperture radii. 
    for ap in ap_radii:
        print('Doing aperture photometry for {}, aperture radius = {} pix, inner annulus radius = {} pix, outer annulus radius = {} pix.'.format(target, ap, an_in, an_out))

        #Declare a new dataframe to hold the information for all targets for this aperture. 
        columns = ['Filename', 'Time UT', 'Time JD', 'Airmass', 'Seeing']
        for i in range(0, len(source_names)):
            columns.append(source_names[i]+' Flux')
            columns.append(source_names[i]+' Flux Error')

        ap_df = pd.DataFrame(index=range(len(reduced_files)), columns=columns)
        output_filename = pines_path/('Objects/'+short_name+'/aper_phot/'+short_name+'_aper_phot_'+str(float(ap))+'_pix_radius.csv')

        #Loop over all images.
        for j in range(len(reduced_files)):
            data = fits.open(reduced_files[j])[0].data
            
            #Read in some supporting information.
            log_path = pines_path/('Logs/'+reduced_files[j].name.split('.')[0]+'_log.txt')
            log = pines_log_reader(log_path)
            header = fits.open(reduced_files[j])[0].header
            date_obs = header['DATE-OBS']
            #Catch a case that can cause datetime strptime to crash; Mimir headers sometimes have DATE-OBS with seconds specified as 010.xx seconds, when it should be 10.xx seconds. 
            if len(date_obs.split(':')[-1].split('.')[0]) == 3:
                date_obs = date_obs.split(':')[0] + ':' + date_obs.split(':')[1] + ':' + date_obs.split(':')[-1][1:]
            
            #Keep a try/except clause here in case other unknown DATE-OBS formats pop up. 
            try:
                date = datetime.datetime.strptime(date_obs, '%Y-%m-%dT%H:%M:%S.%f')
            except:
                print('Header DATE-OBS format does not match the format code in strptime! Inspect/correct the DATE-OBS value.')
                pdb.set_trace()
            days = date.day + hmsm_to_days(date.hour,date.minute,date.second,date.microsecond)
            jd = date_to_jd(date.year,date.month,days)
            ap_df['Filename'][j] = reduced_files[j].name
            ap_df['Time UT'][j] = header['DATE-OBS']
            ap_df['Time JD'][j] = jd
            ap_df['Airmass'][j] = header['AIRMASS']
            ap_df['Seeing'][j] = log['X seeing'][np.where(log['Filename'] == reduced_files[j].name.split('_')[0]+'.fits')[0][0]]
            
            #Get the source positions in this image.
            positions = []
            for i in range(len(source_names)):
                positions.append((centroided_sources[source_names[i]+' X'][j], centroided_sources[source_names[i]+' Y'][j]))

            #Create an aperture centered on this position with radius = ap. 
            apertures = CircularAperture(positions, r=ap)

            #Create an annulus centered on this position. 
            annuli = CircularAnnulus(positions, r_in=an_in, r_out=an_out)

            photometry_tbl = iraf_style_photometry(apertures, annuli, (data*8.21))

            for i in range(len(photometry_tbl)):
                ap_df[source_names[i]+' Flux'][j] = photometry_tbl['flux'][i] 
                ap_df[source_names[i]+' Flux Error'][j] = photometry_tbl['flux_error'][i]
           
           
            # #Get the raw flux in the aperture and annulus. 
            # raw_ap_flux = aperture_photometry(data, apertures)

            #stats = sigma_clipped_stats(data)
            #plt.imshow(8.21*(data - stats[1])/60, origin='lower', vmin=0, vmax=7*8.21*stats[2]/60)

            # annuli_masks = annuli.to_mask(method='center')
            # for i in range(len(annuli_masks)):
            #     annulus_data = annuli_masks[i].multiply(data)
            #     mask = annuli_masks[i].data
            #     annulus_data_1d = annulus_data[mask > 0]
            #     stats = sigma_clipped_stats(annulus_data_1d)
            #     background = stats[1] * apertures[i].area
            #     ap_counts = raw_ap_flux['aperture_sum'][i] - background
            #     ap_df[source_names[i]+' Aper Phot'][j] = ap_counts
            #     ap_df[source_names[i]+' Background'][j] = stats[1]

        #plt.errorbar(ap_df['Time JD'],ap_df[source_names[0]+' Flux']/np.median(ap_df[source_names[0]+' Flux']),yerr=ap_df[source_names[0]+' Flux Error']/np.median(ap_df[source_names[i]+' Flux']), marker='.')
        print('Saving ap = {} aperture photometry output to {}.'.format(ap,output_filename))
        with open(output_filename, 'w') as f:
            for j in range(len(ap_df)):
                if j == 0:
                    f.write('{:>21s}, {:>22s}, {:>17s}, {:>7s}, {:>7s}, '.format('Filename', 'Time UT', 'Time JD', 'Airmass', 'Seeing'))
                    for i in range(len(source_names)):
                        if i != len(source_names) - 1:
                            f.write('{:>20s}, {:>26s}, '.format(source_names[i]+' Flux', source_names[i]+' Flux Error'))
                        else:
                            f.write('{:>20s}, {:>26s}\n'.format(source_names[i]+' Flux', source_names[i]+' Flux Error'))

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
                for i in range(len(source_names)):                    
                    if i != len(source_names) - 1:
                        format_string = '{:20.11f}, {:26.11f}, '
                    else:
                        format_string = '{:20.11f}, {:26.11f}\n'
                    
                    f.write(format_string.format(ap_df[source_names[i]+' Flux'][j], ap_df[source_names[i]+' Flux Error'][j]))
    print('')
    return 