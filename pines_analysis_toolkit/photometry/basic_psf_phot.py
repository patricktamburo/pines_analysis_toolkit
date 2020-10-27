import pines_analysis_toolkit as pat
import pdb
import matplotlib.pyplot as plt
import numpy as np
import natsort
from pathlib import Path
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable
from photutils.psf import BasicPSFPhotometry, IntegratedGaussianPRF, DAOGroup, extract_stars
from astropy.table import Table
from astropy.stats import gaussian_sigma_to_fwhm
from photutils.background import MMMBackground
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import sigma_clipped_stats
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import animation
from astropy.visualization import simple_norm
from photutils import EPSFBuilder, aperture_photometry, EPSFFitter
from astropy.nddata import NDData
from scipy.stats import sigmaclip
from photutils import CircularAperture, CircularAnnulus
import pandas as pd
import datetime
import math 
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
import shutil
import os
import warnings
from astropy.utils.exceptions import AstropyUserWarning
from pines_analysis_toolkit.utils.gif_utils import gif_maker
from photutils.detection import IRAFStarFinder
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from photutils.background import MMMBackground, MADStdBackgroundRMS
from astropy.modeling.fitting import LevMarLSQFitter
from photutils.psf import IterativelySubtractedPSFPhotometry
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.modeling import models, fitting
from scipy import optimize

warnings.simplefilter("ignore", category=AstropyUserWarning)

def basic_psf_phot(target, centroided_sources, plots=False):
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

    def gaussian(p, x, y):
        height, center_x, center_y, width_x, width_y, rotation = p
        rotation = np.deg2rad(rotation)
        x_0 = center_x * np.cos(rotation) - center_y * np.sin(rotation)
        y_0 = center_x * np.sin(rotation) + center_y * np.cos(rotation)

        def rotgauss(x,y):
            xp = x * np.cos(rotation) - y * np.sin(rotation)
            yp = x * np.sin(rotation) + y * np.cos(rotation)
            g = height*np.exp(
                -(((x_0-xp)/width_x)**2+
                  ((y_0-yp)/width_y)**2)/2.)
            return g
        
        g = rotgauss(x,y)

        return g

    def moments(data):
        total = np.nansum(data)
        X, Y = np.indices(data.shape)
        center_x = int(np.shape(data)[1]/2)
        center_y = int(np.shape(data)[0]/2)
        row = data[int(center_x), :]
        col = data[:, int(center_y)]
        width_x = np.nansum(np.sqrt(abs((np.arange(col.size)-center_y)**2*col))
                            /np.nansum(col))
        width_y = np.nansum(np.sqrt(abs((np.arange(row.size)-center_x)**2*row))
                            /np.nansum(row))
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
    short_name = short_name_creator(target)
    reduced_path = pines_path/('Objects/'+short_name+'/reduced/')
    reduced_files = np.array(natsort.natsorted([x for x in reduced_path.glob('*.fits')]))

    centroided_sources.columns = centroided_sources.columns.str.strip()
    source_names = natsort.natsorted(list(set([i[0:-2].replace('X','').replace('Y','').rstrip().lstrip() for i in centroided_sources.keys()])))    

    #Declare a new dataframe to hold the information for all targets for this .
    columns = ['Filename', 'Time UT', 'Time JD', 'Airmass', 'Seeing']
    for i in range(0, len(source_names)):
        columns.append(source_names[i]+' Flux')
        columns.append(source_names[i]+' Flux Error')
    psf_df = pd.DataFrame(index=range(len(reduced_files)), columns=columns)
    output_filename = pines_path/('Objects/'+short_name+'/psf_phot/'+short_name+'_psf_phot.csv')

    for i in range(len(reduced_files)):
        #Read in image data/header. 
        file = reduced_files[i]
        data = fits.open(file)[0].data
        header = fits.open(file)[0].header
        print('{}, image {} of {}.'.format(file.name, i+1, len(reduced_files)))

        #Read in some supporting information.
        log_path = pines_path/('Logs/'+file.name.split('.')[0]+'_log.txt')
        log = pines_log_reader(log_path)
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
        psf_df['Filename'][i] = file.name
        psf_df['Time UT'][i] = header['DATE-OBS']
        psf_df['Time JD'][i] = jd
        psf_df['Airmass'][i] = header['AIRMASS']
        psf_df['Seeing'][i] = log['X seeing'][np.where(log['Filename'] == file.name.split('_')[0]+'.fits')[0][0]]
        
        #Read in source centroids for this image
        x = np.zeros(len(source_names))
        y = np.zeros(len(source_names))
        seeing = psf_df['Seeing'][i]

        for j in range(len(source_names)):
            source = source_names[j]
            x[j] = centroided_sources[source+' X'][i]
            y[j] = centroided_sources[source+' Y'][i]

         #The extract_stars() function requires the input data as an NDData object. 
        nddata = NDData(data=data)  

        #Create table of good star positions
        stars_tbl = Table()
        stars_tbl['x'] = x
        stars_tbl['y'] = y

        size = 25
        x, y = np.meshgrid(np.arange(0,size), np.arange(0,size))

        #Extract star cutouts.
        stars = extract_stars(nddata, stars_tbl, size=size)  

        fitter = fitting.LevMarLSQFitter()

        fig, ax = plt.subplots(nrows=len(stars), ncols=3, sharex=True, sharey=True, figsize=(12,40))

        #Fit a 2D Gaussian to each star. 
        for j in range(len(stars)): 
            star = stars[j]
            source = source_names[j]
            mmm_bkg = MMMBackground()
            cutout = star.data - mmm_bkg(star.data)            

            #Get the star's centroid position in the cutout. 
            dtype = [('x_0', 'f8'), ('y_0', 'f8')]
            pos = Table(data=np.zeros(1, dtype=dtype))
            source_x = stars_tbl['x'][j]
            source_y = stars_tbl['y'][j]
            pos['x_0'] = source_x - int(source_x - size/2 + 1)
            pos['y_0'] = source_y - int(source_y - size/2 + 1)

            parameters = fitgaussian(cutout)
            g2d_fit = gaussian(parameters, x, y)

            avg, med, std = sigma_clipped_stats(cutout)
            im = ax[j,0].imshow(cutout, origin='lower', vmin=med-std, vmax=med+8*std)
            divider = make_axes_locatable(ax[j,0])
            cax = divider.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(im, cax=cax, orientation='vertical')
            ax[j,0].plot(pos['x_0'], pos['y_0'], 'rx')
            ax[j,0].set_ylabel(source)
            ax[j,0].text(pos['x_0'], pos['y_0']+1, '('+str(np.round(source_x,1))+', '+str(np.round(source_y,1))+')', color='r', ha='center')
            ax[j,0].axis('off')

            axins = ax[j,0].inset_axes([0.75, 0.75, 0.25, 0.25])
            axins.set_yticklabels([])
            axins.set_yticks([])
            axins.set_xticklabels([])
            axins.set_xticks([])
            axins.imshow(data, origin='lower', vmin=med-std, vmax=med+8*std)
            axins.plot(source_x, source_y, 'rx')

            im = ax[j,1].imshow(g2d_fit, origin='lower', vmin=med-std, vmax=med+8*std)
            divider = make_axes_locatable(ax[j,1])
            cax = divider.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(im, cax=cax, orientation='vertical')
            ax[j,1].axis('off')

            avg, med, std = sigma_clipped_stats(cutout - g2d_fit)
            im = ax[j,2].imshow(cutout - g2d_fit, origin='lower', vmin=med-std, vmax=med+8*std)
            divider = make_axes_locatable(ax[j,2])
            cax = divider.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(im, cax=cax, orientation='vertical')
            ax[j,2].axis('off')

            if j == 0:
                ax[j,0].set_title('Data')
                ax[j,1].set_title('2D Gaussian Model')
                ax[j,2].set_title('Data - Model')

            plt.tight_layout()

        output_filename = pines_path/('Objects/'+short_name+'/basic_psf_phot/'+reduced_files[i].name.split('_')[0]+'_'+'source_modeling.pdf')
        plt.savefig(output_filename)
        plt.close()

        
    return
           
