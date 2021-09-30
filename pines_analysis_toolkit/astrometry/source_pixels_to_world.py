import pandas as pd 
import pines_analysis_toolkit as pat 
from astropy.io import fits
from astropy.wcs import WCS
import pdb 
import numpy as np 

def source_pixels_to_world(target, source_detect_image_path, force_output_path=''): 
    '''
        Authors:
            Pat Tamburo, Boston University, May 2021
        Purpose: 
            Gets world coordinates of tracked sources (target + references) in the source_detect_image. 
        Inputs:
            target (str): The target's long 2MASS name (e.g., '2MASS J09161504+2139512')
            source_detect_image_path (pathlib.Path): Path to the source detection image used to generate the target_and_references_source_detection.csv file
        Outputs:
            sources (pd.DataFrame): Data frame containing source names, x/y pixel positions, and RA/Decs in the source_detect_image.
            force_output_path (path): if you want to manually set an output directory, specify the top-level here (i.e. the directory containing Logs, Objects, etc.)
    '''

    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pat.utils.pines_dir_check()

    short_name = pat.utils.short_name_creator(target)
    sources_csv_path = pines_path/('Objects/'+short_name+'/sources/target_and_references_source_detection.csv')
    sources_df = pd.read_csv(sources_csv_path)

    #Get the WCS information. 
    source_detect_image = fits.open(source_detect_image_path)
    w = WCS(source_detect_image[0].header)

    #Get the world coordinates of all sources
    source_ras = np.zeros(len(sources_df))
    source_decs = np.zeros(len(sources_df))
    for i in range(len(sources_df)):
        source_x = sources_df['Source Detect X'][i]
        source_y = sources_df['Source Detect Y'][i]
        sky = w.pixel_to_world(source_x, source_y)
        source_ras[i] = sky.ra.value
        source_decs[i] = sky.dec.value
    
    #Add world coordinates to source detection df. 
    sources_df['Source Detect RA'] = source_ras
    sources_df['Source Detect Dec'] = source_decs

    #Write out to update the source detection csv. 
    sources_df.to_csv(sources_csv_path, index=0)
    return sources_df