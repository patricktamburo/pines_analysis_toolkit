from photutils.centroids import centroid_2dg
from photutils.centroids import centroid_1dg
from photutils.centroids import centroid_com
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader
from pines_analysis_toolkit.utils.quick_plot import quick_plot
import natsort
from astropy.stats import sigma_clipped_stats
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pdb
import pandas as pd
from scipy.stats import sigmaclip
import time 
import shutil
import os
import warnings
from astropy.convolution import interpolate_replace_nans
from astropy.utils.exceptions import AstropyUserWarning
from astropy.convolution import Gaussian2DKernel
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

def centroider(target, sources, plots=False, restore=False):
    t1 = time.time()
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    plt.ion()

    kernel = Gaussian2DKernel(x_stddev=1) #For fixing nans in cutouts.

    #If restore == True, read in existing output and return. 
    if restore: 
        centroid_df = pd.read_csv(pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv'),converters={'X Centroids':eval, 'Y Centroids':eval})
        print('Restoring centroider output from {}.'.format(pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv')))
        print('')
        return centroid_df

    #Create subdirectories in sources folder to contain output plots. 
    if plots:
        for name in sources['Name']:
            #If the folders are already there, delete them. 
            source_path = (pines_path/('Objects/'+short_name+'/sources/'+name+'/'))
            if source_path.exists():
                shutil.rmtree(source_path)
            #Create folders.
            os.mkdir(source_path)

    np.seterr(divide='ignore', invalid='ignore') #Suppress some warnings we don't care about in median combining. 

    #Get list of reduced files for target. 
    reduced_path = pines_path/('Objects/'+short_name+'/reduced')
    reduced_files = np.array(natsort.natsorted([x for x in reduced_path.glob('*.fits')]))

    #Declare a new dataframe to hold the centroid information for all sources we want to track. 
    columns = []
    for i in range(0, len(sources)):
        columns.append(sources['Name'][i]+' X')
        columns.append(sources['Name'][i]+' Y')

    centroid_df = pd.DataFrame(index=range(len(reduced_files)), columns=columns)

    log_path = pines_path/('Logs/')
    log_dates = np.array(natsort.natsorted([x.name.split('_')[0] for x in log_path.glob('*.txt')]))

    #Make sure we have logs for all the nights of these data. Need them to account for image shifts. 
    nights = list(set([i.name.split('.')[0] for i in reduced_files]))
    for i in nights:
        if i not in log_dates:
            print('ERROR: {} not in {}. Download it from the PINES server.'.format(i+'_log.txt',log_path))
            pdb.set_trace()

    box_w = 15 #1/2 wdth of the box to cut out around the target position. 
    shift_tolerance = 3 #Number of pixels that the measured centroid can be away from the expected position before trying other centroiding algorithms. 
    for i in range(len(sources)):
        #Get the initial source position.
        x_pos = sources['Source Detect X'][i] 
        y_pos = sources['Source Detect Y'][i] 
        print('')
        print('Getting centroids for {}, source {} of {}.'.format(sources['Name'][i],i+1,len(sources)))
        print('')
        for j in range(len(reduced_files)):
            file = reduced_files[j]
            if i == 0:
                print('{}, {}, image {} of {}.'.format(sources['Name'][i], file.name, j+1, len(reduced_files)))
            else:
                print('{} of {}, {}, image {} of {}.'.format(sources['Name'][i], len(sources)-1, file.name, j+1, len(reduced_files)))

            image = fits.open(file)[0].data

            #Get the measured image shift for this image. 
            log = pines_log_reader(log_path/(file.name.split('.')[0]+'_log.txt'))
            log_ind = np.where(log['Filename'] == file.name.split('_')[0]+'.fits' )[0][0]
            x_shift = float(log['X shift'][log_ind])
            y_shift = float(log['Y shift'][log_ind])

            if np.isnan(x_shift) or np.isnan(y_shift):
                x_shift = 0
                y_shift = 0
            
            #If there are clouds, shifts could have been erroneously high...just zero them?
            if abs(x_shift) > 30:
                x_shift = 0
            if abs(y_shift) > 30:
                y_shift = 0

            #Apply the shift. 
            x_pos = sources['Source Detect X'][i] - x_shift
            y_pos = sources['Source Detect Y'][i] + y_shift

            #Get sigma_clipped_stats of the box around this guess position.
            #stats = sigma_clipped_stats(image[int(y_pos-box_w):int(y_pos+box_w), int(x_pos-box_w):int(x_pos+box_w)])
            cutout = interpolate_replace_nans(image[int(y_pos-box_w):int(y_pos+box_w), int(x_pos-box_w):int(x_pos+box_w)], kernel)

            vals, lower, upper = sigmaclip(cutout, low=1.5, high=2.5) #~3x faster than astropy's sigma_clipped_stats
            med = np.nanmedian(vals)
            std = np.nanstd(vals)

            # #Mask out the image except for around the source position. 
            # mask = np.ones(image.shape, dtype=bool)
            # mask[int(y_pos-box_w):int(y_pos+box_w), int(x_pos-box_w):int(x_pos+box_w)] = False

            # #Do a sigma clipping elimination of pixels around the target, and mask them out. 
            # other_bad_y = np.where(cutout < lower)[1] + int(y_pos) - box_w
            # other_bad_x = np.where(cutout < lower)[0] + int(x_pos) - box_w
            # for k in range(len(other_bad_x)):
            #     mask[other_bad_y[k],[other_bad_x[k]]] = True
            # if len(other_bad_x) == (2*box_w)**2:
            #     plt.imshow(mask, origin='lower')
            #     pdb.set_trace()
            
            
            #Get the centroid. Default to centroid_com, it's the fastest.
            #centroid_x, centroid_y = centroid_com(image-med, mask=mask)

            #If the centroid position exceeds 5 pixels of the x/y pixel position, try a 1d Gaussian detection
            #if (abs(centroid_x - x_pos) > 5) or (abs(centroid_y - y_pos) > 5):
            #    print('centroid_com returned large centroid deviation, trying centroid_1dg')
            
            centroid_x, centroid_y = centroid_2dg(cutout - med, error=1/(np.sqrt(cutout)**2))
            centroid_x += int(x_pos) - box_w
            centroid_y += int(y_pos) - box_w

            # fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10,5), sharex=True, sharey=True)
            # ax[0].imshow(mask, origin='lower')
            # ax[0].set_xlim(int(x_pos-box_w), int(x_pos+box_w))
            # ax[0].set_ylim(int(y_pos-box_w), int(y_pos+box_w))
            # ax[1].imshow(image, origin='lower', vmin=lower, vmax=upper)
            # ax[1].plot(centroid_x, centroid_y, 'rx')

            # pdb.set_trace()
            # plt.close()

            #If that fails, try a 2d Gaussian detection. 
            if (abs(centroid_x - x_pos) > shift_tolerance) or (abs(centroid_y - y_pos) > shift_tolerance):
                #Try without error weighting. 
                centroid_x, centroid_y = centroid_1dg(image)
                if (abs(centroid_x - x_pos) > shift_tolerance) or (abs(centroid_y - y_pos) > shift_tolerance):
                    centroid_x, centroid_y = centroid_2dg(image-med, error=1/(np.sqrt(image)**2))
                    if (abs(centroid_x - x_pos) > shift_tolerance) or (abs(centroid_y - y_pos) > shift_tolerance):
                        centroid_x, centroid_y = centroid_2dg(image)
                        if (abs(centroid_x - x_pos) > shift_tolerance) or (abs(centroid_y - y_pos) > shift_tolerance):
                            centroid_x, centroid_y = centroid_com(image - med)
                            if (abs(centroid_x - x_pos) > shift_tolerance) or (abs(centroid_y - y_pos) > shift_tolerance):
                                centroid_x, centroid_y = centroid_com(image)
                                if (abs(centroid_x - x_pos) > shift_tolerance) or (abs(centroid_y - y_pos) > shift_tolerance):
                                    print('WARNING: large centroid deviation measured, returning predicted position')
                                    print('')
                                    centroid_x = x_pos
                                    centroid_y = y_pos
            #Check that your measured position is actually on the detector. 
            if (centroid_x < 0) or (centroid_y < 0) or (centroid_x > 1023) or (centroid_y > 1023):
                pdb.set_trace()
            
            #Check to make sure you didn't measure nan's. 
            if np.isnan(centroid_x):
                centroid_x = x_pos
                print('NaN returned from centroid algorithm, defaulting to target position in source_detct_image.')
            if np.isnan(centroid_y):
                centroid_y = y_pos
                print('NaN returned from centroid algorithm, defaulting to target position in source_detct_image.')
            

            

            #Record the position.
            centroid_df[sources['Name'][i]+' X'][j] = centroid_x
            centroid_df[sources['Name'][i]+' Y'][j] = centroid_y

            if plots: 
                #Plot
                lock_x = int(centroid_df[sources['Name'][i]+' X'][0])
                lock_y = int(centroid_df[sources['Name'][i]+' Y'][0])
                plt.imshow(image, origin='lower', vmin=med, vmax=med+5*std)
                plt.plot(centroid_x, centroid_y, 'rx')
                plt.ylim(lock_y-box_w,lock_y+box_w-1)
                plt.xlim(lock_x-box_w,lock_x+box_w-1)
                plt.title('CENTROID DIAGNOSTIC PLOT\n'+sources['Name'][i]+', '+reduced_files[j].name+' (image '+str(j+1)+' of '+str(len(reduced_files))+')')
                plt.text(centroid_x, centroid_y+0.5, '('+str(np.round(centroid_x,1))+', '+str(np.round(centroid_y,1))+')', color='r', ha='center')
                plot_output_path = (pines_path/('Objects/'+short_name+'/sources/'+sources['Name'][i]+'/'+str(j).zfill(4)+'.jpg'))
                plt.savefig(plot_output_path)
                plt.close()

    output_filename = pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv')
    #centroid_df.to_csv(pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv'))
    
    print('Saving centroiding output to {}.'.format(output_filename))
    with open(output_filename, 'w') as f:
        for j in range(len(centroid_df)):
            if j == 0:
                for i in range(len(sources['Name'])):
                    if i != len(sources['Name']) - 1:
                        f.write('{:<17s}, {:<17s}, '.format(sources['Name'][i]+' X', sources['Name'][i]+' Y'))
                    else:
                        f.write('{:<17s}, {:<17s}\n'.format(sources['Name'][i]+' X', sources['Name'][i]+' Y'))

            for i in range(len(sources['Name'])):                    
                if i != len(sources['Name']) - 1:
                    format_string = '{:<17.4f}, {:<17.4f}, '
                else:
                    format_string = '{:<17.4f}, {:<17.4f}\n'
                
                f.write(format_string.format(centroid_df[sources['Name'][i]+' X'][j], centroid_df[sources['Name'][i]+' Y'][j]))
    np.seterr(divide='warn', invalid='warn') 
    print('')
    print('centroider runtime: {:.2f} minutes.'.format((time.time()-t1)/60))
    print('')
    return centroid_df