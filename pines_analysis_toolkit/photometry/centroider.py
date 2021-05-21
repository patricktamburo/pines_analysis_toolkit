from photutils.centroids import centroid_2dg
from photutils.centroids import centroid_1dg
from photutils.centroids import centroid_com
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader
from pines_analysis_toolkit.utils.quick_plot import quick_plot
from pines_analysis_toolkit.utils.gif_utils import gif_maker
import natsort
from astropy.stats import sigma_clipped_stats
import numpy as np
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
import pdb
import pandas as pd
from scipy.stats import sigmaclip
import time 
import shutil
import os
import warnings
from astropy.utils.exceptions import AstropyUserWarning
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from astropy.visualization import ZScaleInterval, ImageNormalize, SquaredStretch
from progressbar import ProgressBar
from glob import glob 
from photutils import make_source_mask, CircularAperture
from pines_analysis_toolkit.utils.quick_plot import quick_plot as qp 
import datetime 
import julian 

#Turn off warnings from Astropy because they're incredibly annoying. 
warnings.simplefilter("ignore", category=AstropyUserWarning)

'''Authors:
		Patrick Tamburo, Boston University, June 2020
	Purpose:
        Measures positions of sources in a set of reduced images.
	Inputs:
        target (str): The target's full 2MASS name
        sources (pandas dataframe): List of source names, x and y positions. 
        plots (bool, optional): Whether or not to output plots showing centroid positions. Images output to sources directory within the object directory. 
        restore (bool, optional): Whether or not to restore centroider output that already exists. 
    Outputs:
		centroid_df (pandas DataFrame): X and Y centroid positions for each source. 
	TODO:
        Grab logs automatically? 
        Flag bad centroids?
'''

def centroider(target, sources, output_plots=False, gif=False, restore=False, box_w=8):
    matplotlib.use('TkAgg')
    plt.ioff()
    t1 = time.time()
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)

    kernel = Gaussian2DKernel(x_stddev=1) #For fixing nans in cutouts.

    #If restore == True, read in existing output and return. 
    if restore: 
        centroid_df = pd.read_csv(pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv'),converters={'X Centroids':eval, 'Y Centroids':eval})
        print('Restoring centroider output from {}.'.format(pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv')))
        print('')
        return centroid_df

    #Create subdirectories in sources folder to contain output plots. 
    if output_plots:
        subdirs = glob(str(pines_path/('Objects/'+short_name+'/sources'))+'/*/')
        #Delete any source directories that are already there. 
        for name in subdirs:
            shutil.rmtree(name)

        #Create new source directories.
        for name in sources['Name']:
            source_path = (pines_path/('Objects/'+short_name+'/sources/'+name+'/'))
            os.mkdir(source_path)

    #Read in extra shifts, in case the master image wasn't used for source detection.
    extra_shift_path = pines_path/('Objects/'+short_name+'/sources/extra_shifts.txt')
    extra_shifts = pd.read_csv(extra_shift_path, delimiter=' ', names=['Extra X shift', 'Extra Y shift'])
    extra_x_shift = extra_shifts['Extra X shift'][0]
    extra_y_shift = extra_shifts['Extra Y shift'][0]

    np.seterr(divide='ignore', invalid='ignore') #Suppress some warnings we don't care about in median combining. 

    #Get list of reduced files for target. 
    reduced_path = pines_path/('Objects/'+short_name+'/reduced')
    reduced_filenames = natsort.natsorted([x.name for x in reduced_path.glob('*red.fits')])
    reduced_files = np.array([reduced_path/i for i in reduced_filenames])
    
    #Declare a new dataframe to hold the centroid information for all sources we want to track. 
    columns = []
    columns.append('Filename')
    columns.append('Seeing')
    columns.append('Time (JD UTC)')
    columns.append('Airmass')

    #Add x/y positions and cenroid flags for every tracked source
    for i in range(0, len(sources)):
        columns.append(sources['Name'][i]+' Image X')
        columns.append(sources['Name'][i]+' Image Y')
        columns.append(sources['Name'][i]+' Cutout X')
        columns.append(sources['Name'][i]+' Cutout Y')
        columns.append(sources['Name'][i]+' Centroid Warning')

    centroid_df = pd.DataFrame(index=range(len(reduced_files)), columns=columns)

    log_path = pines_path/('Logs/')
    log_dates = np.array(natsort.natsorted([x.name.split('_')[0] for x in log_path.glob('*.txt')]))

    #Make sure we have logs for all the nights of these data. Need them to account for image shifts. 
    nights = list(set([i.name.split('.')[0] for i in reduced_files]))
    for i in nights:
        if i not in log_dates:
            print('ERROR: {} not in {}. Download it from the PINES server.'.format(i+'_log.txt',log_path))
            pdb.set_trace()

    shift_tolerance = 2.0 #Number of pixels that the measured centroid can be away from the expected position in either x or y before trying other centroiding algorithms. 
    for i in range(len(sources)):
        #Get the initial source position.
        x_pos = sources['Source Detect X'][i] 
        y_pos = sources['Source Detect Y'][i] 
        print('')
        print('Getting centroids for {}, ({:3.1f}, {:3.1f}) in source detection image. Source {} of {}.'.format(sources['Name'][i], x_pos, y_pos, i+1,len(sources)))
        if output_plots:
            print('Saving centroid plots to {}.'.format(pines_path/('Objects/'+short_name+'/sources/'+sources['Name'][i]+'/')))
        pbar = ProgressBar()
        for j in pbar(range(len(reduced_files))):
            centroid_df[sources['Name'][i]+' Centroid Warning'][j] = 0
            file = reduced_files[j]
            image = fits.open(file)[0].data
            #Get the measured image shift for this image. 
            log = pines_log_reader(log_path/(file.name.split('.')[0]+'_log.txt'))
            log_ind = np.where(log['Filename'] == file.name.split('_')[0]+'.fits' )[0][0]

            

            x_shift = float(log['X shift'][log_ind])
            y_shift = float(log['Y shift'][log_ind])

            #Save the filename for readability. Save the seeing for use in variable aperture photometry. Save the time for diagnostic plots. 
            if i == 0:
                centroid_df['Filename'][j] = file.name.split('_')[0]+'.fits'
                centroid_df['Seeing'][j] = log['X seeing'][log_ind]
                time_str = fits.open(file)[0].header['DATE-OBS']
                
                #Correct some formatting issues that can occur in Mimir time stamps. 
                if time_str.split(':')[-1] == '60.00':
                    time_str = time_str[0:14]+str(int(time_str.split(':')[-2])+1)+':00.00'
                elif time_str.split(':')[-1] == '010.00':
                    time_str = time_str[0:17]+time_str.split(':')[-1][1:]

                centroid_df['Time (JD UTC)'][j] = julian.to_jd(datetime.datetime.strptime(time_str, '%Y-%m-%dT%H:%M:%S.%f'))
                centroid_df['Airmass'][j] = log['Airmass'][log_ind]
            
            nan_flag = False #Flag indicating if you should not trust the log's shifts. Set to true if x_shift/y_shift are 'nan' or > 30 pixels. 

            #If bad shifts were measured for this image, skip. 
            if log['Shift quality flag'][log_ind] == 1:
                continue

            if np.isnan(x_shift) or np.isnan(y_shift):
                x_shift = 0
                y_shift = 0
                nan_flag = True

            #If there are clouds, shifts could have been erroneously high...just zero them?
            if abs(x_shift) > 200:
                #x_shift = 0
                nan_flag = True
            if abs(y_shift) > 200:
                #y_shift = 0
                nan_flag = True

            #Apply the shift. NOTE: This relies on having accurate x_shift and y_shift values from the log. 
            #If they're incorrect, the cutout will not be in the right place. 
            #x_pos = sources['Source Detect X'][i] - x_shift + extra_x_shift
            #y_pos = sources['Source Detect Y'][i] + y_shift - extra_y_shift

            x_pos = sources['Source Detect X'][i] - (x_shift - extra_x_shift)
            y_pos = sources['Source Detect Y'][i] + (y_shift - extra_y_shift)

            #TODO: Make all this its own function. 

            #Cutout around the expected position and interpolate over any NaNs (which screw up source detection).
            cutout = interpolate_replace_nans(image[int(y_pos-box_w):int(y_pos+box_w)+1, int(x_pos-box_w):int(x_pos+box_w)+1], kernel=Gaussian2DKernel(x_stddev=0.5))

            #interpolate_replace_nans struggles with edge pixels, so shave off edge_shave pixels in each direction of the cutout. 
            edge_shave = 1      
            cutout = cutout[edge_shave:len(cutout)-edge_shave,edge_shave:len(cutout)-edge_shave]
            
            vals, lower, upper = sigmaclip(cutout, low=1.5, high=2.5) #Get sigma clipped stats on the cutout
            med = np.nanmedian(vals)
            std = np.nanstd(vals)

            try:
                centroid_x_cutout, centroid_y_cutout = centroid_2dg(cutout - med) #Perform centroid detection on the cutout.
            except:
                pdb.set_trace()

            centroid_x = centroid_x_cutout + int(x_pos) - box_w + edge_shave #Translate the detected centroid from the cutout coordinates back to the full-frame coordinates. 
            centroid_y = centroid_y_cutout + int(y_pos) - box_w + edge_shave

            # if i == 0:
            #     qp(cutout)
            #     plt.plot(centroid_x_cutout, centroid_y_cutout, 'rx')

            #     # qp(image)
            #     # plt.plot(centroid_x, centroid_y, 'rx')
            #     pdb.set_trace()

            #If the shifts in the log are not 'nan' or > 200 pixels, check if the measured shifts are within shift_tolerance pixels of the expected position.
            #   If they aren't, try alternate centroiding methods to try and find it.

            #Otherwise, use the shifts as measured with centroid_1dg. PINES_watchdog likely failed while observing, and we don't expect the centroids measured here to actually be at the expected position.
            if not nan_flag:
                #Try a 2D Gaussian detection.
                if (abs(centroid_x - x_pos) > shift_tolerance) or (abs(centroid_y - y_pos) > shift_tolerance):
                    centroid_x_cutout, centroid_y_cutout = centroid_2dg(cutout - med)
                    centroid_x = centroid_x_cutout + int(x_pos) - box_w
                    centroid_y = centroid_y_cutout + int(y_pos) - box_w

                    #If that fails, try a COM detection.
                    if (abs(centroid_x - x_pos) > shift_tolerance) or (abs(centroid_y - y_pos) > shift_tolerance):
                        centroid_x_cutout, centroid_y_cutout = centroid_com(cutout - med)
                        centroid_x = centroid_x_cutout + int(x_pos) - box_w
                        centroid_y = centroid_y_cutout + int(y_pos) - box_w

                        #If that fails, try masking source and interpolate over any bad pixels that aren't in the bad pixel mask, then redo 1D gaussian detection.
                        if (abs(centroid_x - x_pos) > shift_tolerance) or (abs(centroid_y - y_pos) > shift_tolerance):
                            mask = make_source_mask(cutout, nsigma=4, npixels=5, dilate_size=3)   
                            vals, lo, hi = sigmaclip(cutout[~mask])
                            bad_locs =  np.where((mask == False) & ((cutout > hi) | (cutout < lo)))
                            cutout[bad_locs] = np.nan
                            cutout = interpolate_replace_nans(cutout, kernel=Gaussian2DKernel(x_stddev=0.5))
                            
                            centroid_x_cutout, centroid_y_cutout = centroid_1dg(cutout - med)
                            centroid_x = centroid_x_cutout + int(x_pos) - box_w
                            centroid_y = centroid_y_cutout + int(y_pos) - box_w

                            #Try a 2D Gaussian detection on the interpolated cutout
                            if (abs(centroid_x - x_pos) > shift_tolerance) or (abs(centroid_y - y_pos) > shift_tolerance):
                                centroid_x_cutout, centroid_y_cutout = centroid_2dg(cutout - med)
                                centroid_x = centroid_x_cutout + int(x_pos) - box_w
                                centroid_y = centroid_y_cutout + int(y_pos) - box_w

                                #Try a COM on the interpolated cutout.
                                if (abs(centroid_x - x_pos) > shift_tolerance) or (abs(centroid_y - y_pos) > shift_tolerance):
                                    centroid_x_cutout, centroid_y_cutout = centroid_com(cutout)
                                    centroid_x = centroid_x_cutout + int(x_pos) - box_w
                                    centroid_y = centroid_y_cutout + int(y_pos) - box_w

                                    #Last resort: try cutting off the edge of the cutout. Edge pixels can experience poor interpolation, and this sometimes helps. 
                                    if (abs(centroid_x - x_pos) > shift_tolerance) or (abs(centroid_y - y_pos) > shift_tolerance):
                                        cutout = cutout[1:-1, 1:-1]
                                        centroid_x_cutout, centroid_y_cutout = centroid_1dg(cutout - med)
                                        centroid_x = centroid_x_cutout + int(x_pos) - box_w + 1
                                        centroid_y = centroid_y_cutout + int(y_pos) - box_w + 1

                                        #Try with a 2DG
                                        if (abs(centroid_x - x_pos) > shift_tolerance) or (abs(centroid_y - y_pos) > shift_tolerance):
                                            centroid_x_cutout, centroid_y_cutout = centroid_2dg(cutout - med)
                                            centroid_x = centroid_x_cutout + int(x_pos) - box_w + 1
                                            centroid_y = centroid_y_cutout + int(y_pos) - box_w + 1

                                            #If ALL that fails, report the expected position as the centroid. 
                                            if (abs(centroid_x - x_pos) > shift_tolerance) or (abs(centroid_y - y_pos) > shift_tolerance):
                                                print('WARNING: large centroid deviation measured, returning predicted position')
                                                print('')
                                                centroid_df[sources['Name'][i]+' Centroid Warning'][j] = 1
                                                centroid_x = x_pos
                                                centroid_y = y_pos
                                                #pdb.set_trace()

            #Check that your measured position is actually on the detector. 
            if (centroid_x < 0) or (centroid_y < 0) or (centroid_x > 1023) or (centroid_y > 1023):
                #Try a quick mask/interpolation of the cutout. 
                mask = make_source_mask(cutout, nsigma=3, npixels=5, dilate_size=3)   
                vals, lo, hi = sigmaclip(cutout[~mask])
                bad_locs =  np.where((mask == False) & ((cutout > hi) | (cutout < lo)))
                cutout[bad_locs] = np.nan
                cutout = interpolate_replace_nans(cutout, kernel=Gaussian2DKernel(x_stddev=0.5))
                centroid_x, centroid_y = centroid_2dg(cutout - med)
                centroid_x += int(x_pos) - box_w
                centroid_y += int(y_pos) - box_w
                if (centroid_x < 0) or (centroid_y < 0) or (centroid_x > 1023) or (centroid_y > 1023):
                    print('WARNING: large centroid deviation measured, returning predicted position')
                    print('')
                    centroid_df[sources['Name'][i]+' Centroid Warning'][j] = 1
                    centroid_x = x_pos
                    centroid_y = y_pos
                    #pdb.set_trace()
            
            #Check to make sure you didn't measure nan's. 
            if np.isnan(centroid_x):
                centroid_x = x_pos
                print('NaN returned from centroid algorithm, defaulting to target position in source_detct_image.')
            if np.isnan(centroid_y):
                centroid_y = y_pos
                print('NaN returned from centroid algorithm, defaulting to target position in source_detct_image.')
            
            #Record the image and relative cutout positions.
            centroid_df[sources['Name'][i]+' Image X'][j] = centroid_x
            centroid_df[sources['Name'][i]+' Image Y'][j] = centroid_y
            centroid_df[sources['Name'][i]+' Cutout X'][j] = centroid_x_cutout
            centroid_df[sources['Name'][i]+' Cutout Y'][j] = centroid_y_cutout

            if output_plots: 
                #Plot
                lock_x = int(centroid_df[sources['Name'][i]+' Image X'][0])
                lock_y = int(centroid_df[sources['Name'][i]+' Image Y'][0])
                norm = ImageNormalize(data=cutout, interval=ZScaleInterval())
                plt.imshow(image, origin='lower', norm=norm)
                plt.plot(centroid_x, centroid_y, 'rx')
                ap = CircularAperture((centroid_x, centroid_y), r=5)
                ap.plot(lw=2, color='b')
                plt.ylim(lock_y-30,lock_y+30-1)
                plt.xlim(lock_x-30,lock_x+30-1)
                plt.title('CENTROID DIAGNOSTIC PLOT\n'+sources['Name'][i]+', '+reduced_files[j].name+' (image '+str(j+1)+' of '+str(len(reduced_files))+')', fontsize=10)
                plt.text(centroid_x, centroid_y+0.5, '('+str(np.round(centroid_x,1))+', '+str(np.round(centroid_y,1))+')', color='r', ha='center')
                plot_output_path = (pines_path/('Objects/'+short_name+'/sources/'+sources['Name'][i]+'/'+str(j).zfill(4)+'.jpg'))
                plt.gca().set_axis_off()
                plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
                plt.margins(0,0)
                plt.gca().xaxis.set_major_locator(plt.NullLocator())
                plt.gca().yaxis.set_major_locator(plt.NullLocator())
                plt.savefig(plot_output_path, bbox_inches='tight', pad_inches=0, dpi=150)
                plt.close()

        if gif: 
            gif_path = (pines_path/('Objects/'+short_name+'/sources/'+sources['Name'][i]+'/'))
            gif_maker(path=gif_path, fps=10)

    output_filename = pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv')
    #centroid_df.to_csv(pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv'))
    
    print('Saving centroiding output to {}.'.format(output_filename))
    with open(output_filename, 'w') as f:
        for j in range(len(centroid_df)):
            #Write the header line. 
            if j == 0:
                f.write('{:<17s}, '.format('Filename'))
                f.write('{:<15s}, '.format('Time (JD UTC)'))
                f.write('{:<6s}, '.format('Seeing'))
                f.write('{:<7s}, '.format('Airmass'))
                for i in range(len(sources['Name'])):
                    n = sources['Name'][i]
                    if i != len(sources['Name']) - 1:
                        f.write('{:<23s}, {:<23s}, {:<24s}, {:<24s}, {:<34s}, '.format(n+' Image X', n+' Image Y', n+' Cutout X', n+' Cutout Y',  n+' Centroid Warning'))
                    else:
                        f.write('{:<23s}, {:<23s}, {:<24s}, {:<24s}, {:<34s}\n'.format(n+' Image X', n+' Image Y', n+' Cutout X', n+' Cutout Y',  n+' Centroid Warning'))
            
            #Write in the data lines.
            try:
                f.write('{:<17s}, '.format(centroid_df['Filename'][j]))
                f.write('{:<15.7f}, '.format(centroid_df['Time (JD UTC)'][j]))
                f.write('{:<6.1f}, '.format(float(centroid_df['Seeing'][j])))
                f.write('{:<7.2f}, '.format(centroid_df['Airmass'][j]))
            except:
                pdb.set_trace()

            for i in range(len(sources['Name'])):
                n = sources['Name'][i]                    
                if i != len(sources['Name']) - 1:
                    format_string = '{:<23.4f}, {:<23.4f}, {:<24.4f}, {:<24.4f}, {:<34d}, '
                else:
                    format_string = '{:<23.4f}, {:<23.4f}, {:<24.4f}, {:<24.4f}, {:<34d}\n'
                
                f.write(format_string.format(centroid_df[n+' Image X'][j], centroid_df[n+' Image Y'][j], centroid_df[n+' Cutout X'][j], centroid_df[n+' Cutout Y'][j], centroid_df[n+' Centroid Warning'][j]))
    np.seterr(divide='warn', invalid='warn') 
    print('')
    print('centroider runtime: {:.2f} minutes.'.format((time.time()-t1)/60))
    print('')
    return centroid_df