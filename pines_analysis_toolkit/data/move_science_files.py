from astropy.io import fits
import numpy as np
import pdb
import pathlib
import os
import shutil
import time

#Searches date directories for observations of a target, and moves them to a photometry directory.

def move_science_files(target):
    #Get a "short name" for the target for making a directory on PINES hard drive. 
    split_name = target.split('J')[1]
    short_name = '2MASS '+split_name[0:4]+split_name[8]+split_name[9:13]
    top_level_path = '/Users/tamburo/Documents/Data/PINES/'
    data_path = top_level_path + 'Data Staging/'
    target_files = []

    #Find files in data_path
    files = [f for f in os.listdir(data_path) if (os.path.isfile(os.path.join(data_path, f))) and ('.fits' in f)]

    #Loop over files and find where they match target. 
    for file in files:
        file_path = pathlib.Path(str(data_path)+'/'+file)
        object_name = fits.open(file_path)[0].header['OBJECT']
        if object_name == target:
            print(file_path)
            target_files.append(file_path)
    #End loop over files.
    target_files = np.sort(np.array(target_files))

    print('Found ',len(target_files),' files.')
    time.sleep(2)
    target_directory = pathlib.Path(top_level_path+'Objects/'+short_name)
    if not os.path.exists(target_directory):
        print('Making directory for ',short_name,'.')
        os.mkdir(target_directory)
        print('Making reduction/aperture photometry directories.')
        os.mkdir(pathlib.Path(str(target_directory)+'/aper_phot'))
        os.mkdir(pathlib.Path(str(target_directory)+'/domeflat_reduced'))
        os.mkdir(pathlib.Path(str(target_directory)+'/skyflat_reduced'))
        os.mkdir(pathlib.Path(str(target_directory)+'/raw'))
        os.mkdir(pathlib.Path(str(target_directory)+'/analysis'))


        # txt = input('Moving files to '+str(target_directory)+'/raw, type y to continue ')
        # if txt == 'y':
        target_raw_directory = pathlib.Path(str(target_directory)+'/raw')
        for file in target_files:
            shutil.move(str(file),target_raw_directory)
    else:
        # print(target_directory,' already exists, moving files there!')
        # txt = input('Moving files to '+str(target_directory)+'/raw, type y to continue ')
        # if txt == 'y':
        target_raw_directory = pathlib.Path(str(target_directory)+'/raw')
        for file in target_files:
            try:
                shutil.move(str(file),target_raw_directory)
            except:
                os.remove(str(file))
    print('move_science_files complete!')
