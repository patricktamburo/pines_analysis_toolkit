from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import pdb
import glob 
from astropy.stats import sigma_clipped_stats
import pickle
from photutils import CircularAperture
from photutils import DAOStarFinder
import warnings
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
import natsort
import pandas as pd
def movie(long_target_name,ap_rad, dates_to_examine=[],target_number=0,box_size=300,show_centroid=1):
    import time 
    '''Authors:
            Patrick Tamburo, Boston University, Jan 2020, July 2020
        Purpose:
        Inputs:
        Outputs:
        TODO:
    '''
    warnings.filterwarnings('ignore', category=UserWarning, append=True) #This turns of NaN warnings in sigma_clipped_stats, otherwise we'd get a warning every line. 
    
    pines_path = pines_dir_check()
    target_name = short_name_creator(long_target_name)
    object_path = pines_path/('Objects/'+target_name+'/')
    reduced_path = object_path/'reduced/'
    reduced_files = np.array(natsort.natsorted([x for x in reduced_path.glob('*.fits')]))
    centroid_path = object_path/'sources/'

    positions = pd.read_csv(centroid_path/'target_and_references_centroids.csv')
    x_positions = np.array(positions[positions.keys()[2*target_number]])
    y_positions = np.array(positions[positions.keys()[2*target_number+1]])
    initial_position = (int(x_positions[0]),int(y_positions[0]))

    #Get list of files
    all_frame_list = np.array([])
    for i in range(np.size(reduced_files)):
        file_name = reduced_files[i].name
        all_frame_list = np.append(all_frame_list, file_name)
    num_files = len(reduced_files)

    dates = np.array([reduced_files[i].name.split('.')[0] for i in range(len(reduced_files))])

    plt.ion()
    fig, ax = plt.subplots(1,1,figsize=(8,7))
    for i in range(len(reduced_files)):
        if len(dates_to_examine) != 0:
            if reduced_files[i].name.split('.')[0] in dates_to_examine:
                image_path = reduced_files[i]
                title = image_path.name
                image = fits.open(image_path)[0].data
                header = fits.open(image_path)[0].header
                frame = image[initial_position[1]-int(box_size/2):initial_position[1]+int(box_size/2),initial_position[0]-int(box_size/2):initial_position[0]+int(box_size/2)]
                ap = CircularAperture((x_positions[i],y_positions[i]),r=ap_rad)
                avg,med,std = sigma_clipped_stats(image)
                im = ax.imshow(image,origin='lower',vmin=med,vmax=med+5*std)
                ax.set_xlim(initial_position[0]-int(box_size/2),initial_position[0]+int(box_size/2))
                ax.set_ylim(initial_position[1]-int(box_size/2),initial_position[1]+int(box_size/2))
                cb = fig.colorbar(im, orientation='vertical',label='Counts')
                if show_centroid:
                    ap.plot(color='b',lw=2)
                    ax.plot(x_positions[i],y_positions[i],'bx')
                ax.set_title(title)
                plt.pause(0.01)
                ax.cla()
                cb.remove()
        else:
            image_path = reduced_files[i]
            title = image_path.name
            image = fits.open(image_path)[0].data
            header = fits.open(image_path)[0].header
            frame = image[initial_position[1]-int(box_size/2):initial_position[1]+int(box_size/2),initial_position[0]-int(box_size/2):initial_position[0]+int(box_size/2)]
            ap = CircularAperture((x_positions[i],y_positions[i]),r=ap_rad)
            avg,med,std = sigma_clipped_stats(image)
            im = ax.imshow(image,origin='lower',vmin=med,vmax=med+5*std)
            ax.set_xlim(initial_position[0]-int(box_size/2),initial_position[0]+int(box_size/2))
            ax.set_ylim(initial_position[1]-int(box_size/2),initial_position[1]+int(box_size/2))
            cb = fig.colorbar(im, orientation='vertical',label='Counts')
            if show_centroid:
                ap.plot(color='b',lw=2)
                ax.plot(x_positions[i],y_positions[i],'bx')
            ax.set_title(title)
            plt.pause(0.01)
            ax.cla()
            cb.remove()