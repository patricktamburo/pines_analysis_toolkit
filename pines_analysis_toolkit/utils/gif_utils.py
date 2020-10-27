import imageio
from os import listdir
from os.path import isfile, join
import os
from pathlib import Path
from natsort import natsorted
import numpy as np
import pdb 

''' 
    Author: 
        Patrick Tamburo, July 2020
    Purpose: 
        Creates a gif from a directory of images. 
    Inputs:
        path (str): path to file directory.
        ext (str): extension of images in directory (e.g., '.jpg' or '.png')
        fps (int): frames-per-second of gif. 
        title (str): the title of the gif to save. By default, 'animation'
        delete (bool): whether or not to delete the frames after the gif is made. 
    Outputs:
        Saves title.gif to path.
    TODO:
        Compression of final gif. 
        Implement delete functionality. 
'''

def gif_maker(path='.', ext='.jpg', fps=30, title='animation', delete=False):
    images = []
    file_path = Path(path)
    files = np.array(natsorted([x for x in file_path.glob('*'+ext)]))
    print('Reading in files...')
    for i in range(len(files)):
        print('{}, {} of {}'.format(files[i].name, i+1, len(files)))
        images.append(imageio.imread(files[i]))
    
    print('')
    print('Making gif...')
    if not os.path.exists(file_path/('animation/')):
        os.mkdir(file_path/('animation/'))
    imageio.mimsave(file_path/('animation/'+title+'.gif'), images, fps=fps)
    print('')
    print('Saved gif to {}'.format(file_path/('animation/'+title+'.gif')))
    print('gif size: {} Mb'.format(np.round(os.stat(file_path/('animation/'+title+'.gif')).st_size / (1000**2),1)))

