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
from photutils import EPSFBuilder, aperture_photometry
from astropy.nddata import NDData
from scipy.stats import sigmaclip
from photutils import CircularAperture, CircularAnnulus
import pandas as pd
import datetime
import math 

import warnings
from astropy.utils.exceptions import AstropyUserWarning

#I don't need you to tell me there are NaNs in my data every god damn image. 
warnings.simplefilter("ignore", category=AstropyUserWarning)

def psf_phot(target, centroided_sources):
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
    
    def psf_residual_plot(star, result_tab, residual_image, rotate=False):
        fig = plt.figure(figsize=(7,7))
        ax1 = fig.add_subplot(221)
        avg, med, std = sigma_clipped_stats(star)
        im1 = ax1.imshow(star, origin='lower', vmin=0,vmax=med+5*std)
        ax1.plot(result_tab['x_fit'], result_tab['y_fit'],marker='x',color='tab:blue')
        ax1.set_title('Original image', color='tab:blue')
            
        ax2 = fig.add_subplot(222)
        avg, med, std = sigma_clipped_stats(residual_image)
        im2 = ax2.imshow(residual_image, origin='lower', vmin=0,vmax=med+5*std)
        ax2.plot(result_tab['x_fit'], result_tab['y_fit'],marker='x',color='tab:orange')
        ax2.set_title('PSF extraction residual image', color='tab:orange')

        ax3 = fig.add_subplot(223, projection='3d')
        x = np.arange(0,star.data.shape[0])
        y = np.arange(0,star.data.shape[1])
        x_arr, y_arr = np.meshgrid(x,y)

        avg, med, std = sigma_clipped_stats(star.data)
        ax3.plot_wireframe(x_arr, y_arr, star.data)
        ax3.set_xticks([])
        ax3.set_yticks([])
        ax3.set_zlabel('Counts', fontsize=8)
        
        ax4 = fig.add_subplot(224, projection='3d')
        avg, med, std = sigma_clipped_stats(residual_image)
        ax4.plot_wireframe(x_arr, y_arr, residual_image, color='tab:orange')
        ax4.set_zlim(ax3.get_zlim())
        ax4.set_zlabel('Counts', fontsize=8)
        ax4.set_xticks([])
        ax4.set_yticks([])
        
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cb = fig.colorbar(im1, cax=cax, orientation='vertical', label='Counts')

        divider = make_axes_locatable(ax2)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cb = fig.colorbar(im2, cax=cax, orientation='vertical', label='Counts')

        plt.tight_layout()

        # def rotate(angle):
        #     ax3.view_init(elev=15, azim=angle)
        #     ax4.view_init(elev=15, azim=angle)
        
        # if rotate:
        #     rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0,362,2),interval=50)
        #     rot_animation.save('/Users/tamburo/Desktop/rotation.gif', dpi=80, writer='imagemagick') 

    plt.ion()
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    bp_reduced_path = pines_path/('Objects/'+short_name+'/bp_bg_reduced/')
    bp_reduced_files = np.array(natsort.natsorted([x for x in bp_reduced_path.glob('*.fits')]))

    centroided_sources.columns = centroided_sources.columns.str.strip()
    source_names = natsort.natsorted(list(set([i[0:-2].replace('X','').replace('Y','').rstrip().lstrip() for i in centroided_sources.keys()])))
    if len(source_names) != len(centroided_sources.keys())/2:
        print('ERROR: something going wrong grabbing source names. ')
        return

    #Declare a new dataframe to hold the information for all targets for this .
    columns = ['Filename', 'Time UT', 'Time JD', 'Airmass', 'Seeing']
    for i in range(0, len(source_names)):
        columns.append(source_names[i]+' Flux')
        columns.append(source_names[i]+' Flux Error')
    psf_df = pd.DataFrame(index=range(len(bp_reduced_files)), columns=columns)
    output_filename = pines_path/('Objects/'+short_name+'/psf_phot/'+short_name+'_psf_phot.csv')

    for i in range(len(bp_reduced_files)):
        #Read in image data/header. 
        file = bp_reduced_files[i]
        data = fits.open(file)[0].data
        header = fits.open(file)[0].header
        print('{}, image {} of {}.'.format(file.name, i+1, len(bp_reduced_files)))

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

        #Extract 25 x 25 pixel cutouts of our stars, so let’s explicitly exclude stars that are too close to the image boundaries (because they cannot be extracted).
        size = 25
        hsize = (size - 1) / 2
        mask = ((x > hsize) & (x < (data.shape[1] -1 - hsize)) & (y > hsize) & (y < (data.shape[0] -1 - hsize)) & (y > 100) & (y < 923))

        #Create table of good star positions
        stars_tbl = Table()
        stars_tbl['x'] = x[mask]  
        stars_tbl['y'] = y[mask]  
        
        #Subtract background (star cutouts from which we build the ePSF must have background subtracted).
        mean_val, median_val, std_val = sigma_clipped_stats(data, sigma=2.)  
        data -= median_val

        #The extract_stars() function requires the input data as an NDData object. 
        nddata = NDData(data=data)  

        #Extract star cutouts.
        stars = extract_stars(nddata, stars_tbl, size=size)  

        #Correct NaNs in the star cutouts.  
        for j in range(len(stars)):
            for y_pix in range(0, size):
                for x_pix in range(0, size):
                    if np.isnan(stars[j].data[y_pix, x_pix]):
                        x_start = max(0, x_pix - 1)
                        x_end = min(size-1, x_pix + 1)
                        y_start = max(0, y_pix - 1)
                        y_end = min(size-1, y_pix + 1)

                        cutout_pix = []
                        for yy in range(y_start, y_end + 1):
                            for xx in range(x_start, x_end + 1):
                                if not np.isnan(stars[j].data[yy, xx]):
                                    cutout_pix.append(stars[j].data[yy, xx])
                        stars[j].data[y_pix, x_pix] = np.median(cutout_pix)

        # #Correct outliers in the star cutouts.  
        # sigma = 5.
        # box_size = 3
        # h_size = int(np.floor(box_size/2))
        # for j in range(len(stars)):
        #     for y_pix in range(0, size):
        #         for x_pix in range(0, size):
        #             target_pixel_val = stars[j].data[y_pix, x_pix]
                    
        #             x_start = max(0, x_pix - h_size)
        #             x_end = min(size-h_size, x_pix + h_size)
        #             y_start = max(0, y_pix -h_size)
        #             y_end = min(size-h_size, y_pix + h_size)

        #             cutout_pix = []
        #             for yy in range(y_start, y_end + 1):
        #                 for xx in range(x_start, x_end + 1):
        #                     if stars[j].data[yy, xx] != target_pixel_val:
        #                         cutout_pix.append(stars[j].data[yy, xx])
        #             vals, lo, hi = sigmaclip(cutout_pix, low=2., high=2.)

        #             if abs(target_pixel_val - np.mean(vals)) / np.std(vals) > sigma:
        #                 stars[j].data[y_pix, x_pix] = np.median(cutout_pix)
                        

        #Plot. 
        nrows = 5
        ncols = 5
        fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 10), squeeze=True)
        ax = ax.ravel()
        for j in range(len(stars)):
            norm = simple_norm(stars[j], 'log', percent=99.)
            ax[j].imshow(stars[j], norm=norm, origin='lower', cmap='viridis')

        #Construct the ePSF using the star cutouts. 
        epsf_builder = EPSFBuilder(oversampling=2, maxiters=3, progress_bar=False)  
        epsf, fitted_stars = epsf_builder(stars)  

    
        #Plot the ePSF. 
        plt.figure()
        norm = simple_norm(epsf.data, 'log', percent=99.)
        plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
        cb = plt.colorbar()
        plt.tight_layout()
        pdb.set_trace()

        for j in range(len(stars)):
            star = stars[j]
            source_name = source_names[j]
            sigma_psf = 1.85
            #psf_model = IntegratedGaussianPRF(sigma=sigma_psf)
            #psf_model.x_0.fixed = True
            #psf_model.y_0.fixed = True
            #psf_model.sigma.fixed = False

            N = 1
            dtype = [('x_0', 'f8'), ('y_0', 'f8')]
            pos = Table(data=np.zeros(N, dtype=dtype))
            source_x = stars_tbl['x'][j]
            source_y = stars_tbl['y'][j]
            pos['x_0'] = source_x - int(source_x - size/2 - 1)
            pos['y_0'] = source_y - int(source_y - size/2 - 1)

            daogroup = DAOGroup(4.0*sigma_psf*gaussian_sigma_to_fwhm)
            mmm_bkg = MMMBackground()
            photometry = BasicPSFPhotometry(group_maker=daogroup,
                                bkg_estimator=mmm_bkg,
                                psf_model=epsf,
                                fitter=LevMarLSQFitter(),
                                fitshape=(17,17),
                                aperture_radius=4.)
            

            result_tab = photometry(image=star, init_guesses=pos)
            residual_image = photometry.get_residual_image()
            psf_df[source_name+' Flux'][i] = result_tab['flux_fit'][0]
            psf_df[source_name+' Flux Error'][i] = result_tab['flux_unc'][0]
            if j == 0:
                print('(X, Y): ({}, {})'.format(np.round(result_tab['x_fit'][0],1), np.round(result_tab['y_fit'][0],1)))
    
    pdb.set_trace()
           

if __name__ == '__main__':
    target = '2MASS J23255604-0259508'
    sources = pat.photometry.ref_star_chooser(target, restore=True)
    centroided_sources = pat.photometry.centroider(target, sources, restore=True)
    psf_phot(target, centroided_sources)
