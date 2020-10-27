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

warnings.simplefilter("ignore", category=AstropyUserWarning)

def epsf_phot(target, centroided_sources, plots=False):
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

    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    reduced_path = pines_path/('Objects/'+short_name+'/reduced/')
    reduced_filenames = natsort.natsorted([x.name for x in reduced_path.glob('*.fits')])
    reduced_files = np.array([reduced_path/i for i in reduced_filenames])

    centroided_sources.columns = centroided_sources.columns.str.strip()
    source_names = natsort.natsorted(list(set([i[0:-2].replace('X','').replace('Y','').rstrip().lstrip() for i in centroided_sources.keys()])))
    
    #Create output plot directories for each source.
    if plots:
        for name in source_names:
            #If the folders are already there, delete them. 
            source_path = (pines_path/('Objects/'+short_name+'/psf_phot/'+name+'/'))
            if source_path.exists():
                shutil.rmtree(source_path)
            #Create folders.
            os.mkdir(source_path)

    #Declare a new dataframe to hold the information for all targets for this .
    columns = ['Filename', 'Time UT', 'Time JD', 'Airmass', 'Seeing']
    for i in range(0, len(source_names)):
        columns.append(source_names[i]+' Flux')
        columns.append(source_names[i]+' Flux Error')
    psf_df = pd.DataFrame(index=range(len(reduced_files)), columns=columns)
    output_filename = pines_path/('Objects/'+short_name+'/psf_phot/'+short_name+'_psf_phot.csv')

    for i in range(0, len(reduced_files)):
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
        for j in range(len(source_names)):
            source = source_names[j]
            x[j] = centroided_sources[source+' X'][i]
            y[j] = centroided_sources[source+' Y'][i]

        #Extract pixel cutouts of our stars, so let’s explicitly exclude stars that are too close to the image boundaries (because they cannot be extracted).
        size = 13
        hsize = (size - 1) / 2
        #mask = ((x > hsize) & (x < (data.shape[1] -1 - hsize)) & (y > hsize) & (y < (data.shape[0] -1 - hsize)) & (y > 100) & (y < 923))

        #Create table of good star positions
        stars_tbl = Table()
        stars_tbl['x'] = x
        stars_tbl['y'] = y
        
        #Subtract background (star cutouts from which we build the ePSF must have background subtracted).
        mean_val, median_val, std_val = sigma_clipped_stats(data, sigma=2.)  
        data -= median_val
        
        #Replace nans in data using Gaussian. 
        kernel = Gaussian2DKernel(x_stddev=1)
        data = interpolate_replace_nans(data, kernel)

        #The extract_stars() function requires the input data as an NDData object. 
        nddata = NDData(data=data)  

        #Extract star cutouts.
        stars = extract_stars(nddata, stars_tbl, size=size)  
                        

        #Plot. 
        # nrows = 5
        # ncols = 5
        # fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 10), squeeze=True)
        # ax = ax.ravel()
        # for j in range(len(stars)):           
        #     norm = simple_norm(stars[j], 'log', percent=99.)
        #     ax[j].imshow(stars[j].data, norm=norm, origin='lower', cmap='viridis')

        #Construct the ePSF using the star cutouts.
        epsf_fitter = EPSFFitter()
        epsf_builder = EPSFBuilder(maxiters=4, progress_bar=False, fitter=epsf_fitter)   

        try:
            epsf, fitted_stars = epsf_builder(stars)
            output_filename = pines_path/('Objects/'+short_name+'/psf_phot/'+short_name+'_psf_phot.csv')

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
                                    fitshape=(13,13),
                                    aperture_radius=4.)
                

                result_tab = photometry(image=star, init_guesses=pos)
                residual_image = photometry.get_residual_image()
                psf_df[source_name+' Flux'][i] = result_tab['flux_fit'][0]
                psf_df[source_name+' Flux Error'][i] = result_tab['flux_unc'][0]

                if plots:
                    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(12,4))
                    im = ax[0].imshow(star, origin='lower')
                    divider = make_axes_locatable(ax[0])
                    cax = divider.append_axes('right', size='5%', pad=0.05)
                    fig.colorbar(im, cax=cax, orientation='vertical')
                    ax[0].plot(result_tab['x_fit'][0], result_tab['y_fit'][0], 'rx')
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
                    plt.suptitle(source_name+'\n'+reduced_files[i].name+', image '+str(i+1)+' of '+str(len(reduced_files)))
                    plt.subplots_adjust(wspace=0.5, top=0.95, bottom = 0.05)
                    plot_output_name = pines_path/('Objects/'+short_name+'/psf_phot/'+source_name+'/'+str(i).zfill(4)+'.jpg')
                    plt.savefig(plot_output_name)
                    plt.close()
        except:
            print('')
            print('EPSF BUILDER FAILED, SKIPPING IMAGE.')
            print('')
        #Plot the ePSF. 
        # plt.figure()
        # norm = simple_norm(epsf.data, 'log', percent=99.)
        # plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
        # cb = plt.colorbar()
        # plt.tight_layout()   

        

    print('Saving psf photometry output to {}.'.format(output_filename))
    with open(output_filename, 'w') as f:
        for j in range(len(psf_df)):
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
            if type(psf_df['Seeing'][j]) == str:
                psf_df['Seeing'][j] = float(psf_df['Seeing'][j])

            #Do a try/except clause for writeout, in case it breaks in the future. 
            try:
                f.write(format_string.format(psf_df['Filename'][j], psf_df['Time UT'][j], psf_df['Time JD'][j], psf_df['Airmass'][j], psf_df['Seeing'][j]))
            except:
                print('Writeout failed! Inspect quantities you are trying to write out.')
                pdb.set_trace()
            for i in range(len(source_names)):                    
                if i != len(source_names) - 1:
                    format_string = '{:20.11f}, {:26.11f}, '
                else:
                    format_string = '{:20.11f}, {:26.11f}\n'
                
                f.write(format_string.format(psf_df[source_names[i]+' Flux'][j], psf_df[source_names[i]+' Flux Error'][j]))
    print('')    
    return
           
