from photutils.centroids import centroid_2dg
from photutils.centroids import centroid_1dg
from photutils.centroids import centroid_com
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader
import natsort
from astropy.stats import sigma_clipped_stats
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pdb
import pandas as pd

'''Authors:
		Patrick Tamburo, Boston University, June 2020
	Purpose:
        Measures positions of sources in a set of reduced images.
	Inputs:
        target (str): The target's full 2MASS name
        sources (pandas dataframe): List of source names, x and y positions. 
        plots (bool, optional): Whether or not to plot sources and measured centroids as you go along. 
        restore (bool, optional): Whether or not to restore centroider output that already exists. 
    Outputs:
		x_centroids (numpy array): Array of x positions, num_targets x num_images large. 
        y_centroids (numpy array): Array of y positions, num_targets x num_images large. 
	TODO:
        Grab logs automatically? 
        Flag bad centroids?
'''

def centroider(target, sources, plots=False, restore=False):
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)

    #If restore == True, read in existing output and return. 
    if restore: 
        sources = pd.read_csv(pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv'),converters={'X Centroids':eval, 'Y Centroids':eval}).drop(columns=['Unnamed: 0'])
        print('Restoring centroider output from {}.'.format(pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv')))
        print('')
        return sources

    plt.ion() 
    np.seterr(divide='ignore', invalid='ignore') #Suppress some warnings we don't care about in median combining. 

    #Declare some new empty columns in the sources dataframe.
    sources['X Centroids'] = [[] for x in range(len(sources))]
    sources['Y Centroids'] = [[] for x in range(len(sources))]

    #Get list of reduced files for target. 
    reduced_path = pines_path/('Objects/'+short_name+'/reduced')
    reduced_files = np.array(natsort.natsorted([x for x in reduced_path.glob('*.fits')]))

    log_path = pines_path/('Logs/')
    log_dates = np.array(natsort.natsorted([x.name.split('_')[0] for x in log_path.glob('*.txt')]))

    #Make sure we have logs for all the nights of these data. Need them to account for image shifts. 
    nights = list(set([i.name.split('.')[0] for i in reduced_files]))
    for i in nights:
        if i not in log_dates:
            print('ERROR: {} not in {}. Download it from the PINES server.'.format(i+'_log.txt',log_path))
            pdb.set_trace()

    box_w = 10 #1/2 wdth of the box to cut out around the target position. 
    clip_val  = 4 #Sigma value below which pixels are masked. 
    for i in range(len(sources)):
        #Get the initial source position.
        x_pos = sources['X'][i] 
        y_pos = sources['Y'][i] 
        print('')
        print('Getting centroids for {}, source {} of {}.'.format(sources['Name'][i],i+1,len(sources)))
        print('')
        for j in range(len(reduced_files)):
            file = reduced_files[j]
            image = fits.open(file)[0].data

            #Get the measured image shift for this image. 
            log = pines_log_reader(log_path/(file.name.split('.')[0]+'_log.txt'))
            log_ind = np.where(log['Filename'] == file.name.split('_')[0]+'.fits' )[0][0]
            x_shift = float(log['X shift'][log_ind])
            y_shift = float(log['Y shift'][log_ind])

            #Apply the shift. 
            x_pos = sources['X'][i] - x_shift
            y_pos = sources['Y'][i] + y_shift

            #Get sigma_clipped_stats of the box around this guess position.
            stats = sigma_clipped_stats(image[int(y_pos-box_w):int(y_pos+box_w), int(x_pos-box_w):int(x_pos+box_w)])

            #Mask out the image except for around the source position. 
            mask = np.ones(image.shape, dtype=bool)
            mask[int(y_pos-box_w):int(y_pos+box_w), int(x_pos-box_w):int(x_pos+box_w)] = False

            #Do a sigma clipping elimination of pixels around the target, and mask them out. 
            other_bad_y = np.where(image[int(y_pos-box_w):int(y_pos+box_w), int(x_pos-box_w):int(x_pos+box_w)] < stats[1]+clip_val*stats[2])[0] + int(y_pos) - box_w
            other_bad_x = np.where(image[int(y_pos-box_w):int(y_pos+box_w), int(x_pos-box_w):int(x_pos+box_w)] < stats[1]+clip_val*stats[2])[1] + int(x_pos) - box_w
            for k in range(len(other_bad_x)):
                mask[other_bad_y[k],[other_bad_x[k]]] = True

            #Get the centroid. 
            #centroid_x, centroid_y = centroid_1dg(image-stats[1], error=1/(np.sqrt(image)**2), mask=mask)
            #centroid_x, centroid_y = centroid_2dg(image-stats[1], error=1/(np.sqrt(image)**2) mask=mask)
            centroid_x, centroid_y = centroid_com(image-stats[1], mask=mask)

            if plots: 
                #Plot
                plt.imshow(-1*image*(mask-1), origin='lower', vmin=stats[1], vmax=stats[1]+7*stats[2])
                plt.plot(centroid_x, centroid_y, 'rx')
                plt.ylim(y_pos-box_w,y_pos+box_w)
                plt.xlim(x_pos-box_w,x_pos+box_w)
                plt.title(sources['Name'][i]+', image '+str(j+1)+' of '+str(len(reduced_files)))
                plt.pause(0.01)
                plt.clf()

            #Record the position.
            sources['X Centroids'].iloc[i].extend([centroid_x])
            sources['Y Centroids'].iloc[i].extend([centroid_y])
    sources.to_csv(pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv'))
    np.seterr(divide='warn', invalid='warn') 
    return sources