import pdb
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
import numpy as np
import natsort
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
import datetime 
import math

'''Authors:
		Patrick Tamburo, Boston University, June 2020
	Purpose:
        Performs aperture photometry on a set of reduced images given dataframe of source positions. 
	Inputs:
        target (str): The target's full 2MASS name.
        sources (pandas dataframe): List of source names, x and y positions in every image. 
        ap_radii (list of floats): List of aperture radii in pixels for which aperture photometry wil be performed. 
        an_in (float, optional): The inner radius of the annulus used to estimate background, in pixels. 
        an_out (float, optional): The outer radius of the annulus used to estimate background, in pixels. 
    Outputs:
        None. 
	TODO:
'''

def aper_phot(target, sources, ap_radii, an_in=12., an_out=30.):
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
    
    plt.ion()
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)

    #Get list of reduced files for target. 
    reduced_path = pines_path/('Objects/'+short_name+'/reduced')
    reduced_files = np.array(natsort.natsorted([x for x in reduced_path.glob('*.fits')]))

    #Loop over all aperture radii. 
    for ap in ap_radii:
        print('Doing aperture photometry for {}, aperture radius = {} pix, inner annulus radius = {} pix, outer annulus radius = {} pix.'.format(target, ap, an_in, an_out))

        #Declare a new dataframe to hold the information for all targets for this aperture. 
        columns = ['Time UT', 'Time JD', 'Airmass', sources['Name'][0]+' Aper Phot', sources['Name'][0]+' Background']
        for i in range(1, len(sources)):
            columns.append(sources['Name'][i]+' Aper Phot')
            columns.append(sources['Name'][i]+' Background')

        ap_df = pd.DataFrame(index=range(len(reduced_files)), columns=columns)
        output_filename = pines_path/('Objects/'+short_name+'/aper_phot/'+short_name+'_aper_phot_'+str(float(ap))+'_pix_radius.csv')

        #Loop over all images.
        for j in range(len(reduced_files)):
            data = fits.open(reduced_files[j])[0].data

            #Read in some supporting information.
            header = fits.open(reduced_files[j])[0].header
            date = datetime.datetime.strptime(header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
            days = date.day + hmsm_to_days(date.hour,date.minute,date.second,date.microsecond)
            jd = date_to_jd(date.year,date.month,days)
            ap_df['Time UT'][j] = header['DATE-OBS']
            ap_df['Time JD'][j] = jd
            ap_df['Airmass'][j] = header['AIRMASS']

            #Get the source positions in this image.
            positions = []
            for i in range(len(sources)):
                positions.append((sources['X Centroids'][i][j], sources['Y Centroids'][i][j]))
            
            #Create an aperture centered on this position with radius = ap. 
            apertures = CircularAperture(positions, r=ap)

            #Create an annulus centered on this position. 
            annuli = CircularAnnulus(positions, r_in=an_in, r_out=an_out)

            #Get the raw flux in the aperture and annulus. 
            raw_ap_flux = aperture_photometry(data, apertures)

            # stats = sigma_clipped_stats(data)
            # plt.imshow(data, origin='lower', vmin=stats[1], vmax=stats[1]+7*stats[2])
            # apertures.plot(color='white', lw=2)
            # annuli.plot(color='red', lw=2)                
            annuli_masks = annuli.to_mask(method='center')
            for i in range(len(annuli_masks)):
                annulus_data = annuli_masks[i].multiply(data)
                mask = annuli_masks[i].data
                annulus_data_1d = annulus_data[mask > 0]
                stats = sigma_clipped_stats(annulus_data_1d)
                background = stats[1] * apertures[i].area
                ap_counts = raw_ap_flux['aperture_sum'][i] - background
                ap_df[sources['Name'][i]+' Aper Phot'][j] = ap_counts
                ap_df[sources['Name'][i]+' Background'][j] = stats[1]

        # for i in range(len(sources)):
        #     name = sources['Name'][i] + ' Aper Phot'
        #     plt.plot(ap_df['Time JD'], ap_df[name],marker='.')
        ap_df.to_csv(output_filename, index=False)
        print('Saving ap = {} aperture photometry output to {}.'.format(ap,output_filename))
        print('')
    return 