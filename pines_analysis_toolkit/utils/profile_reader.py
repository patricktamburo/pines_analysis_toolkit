import pandas as pd 
import pdb 
import distutils
import os 
from glob import glob
from natsort import natsorted 
import pines_analysis_toolkit as pat 
import numpy as np 

def profile_reader(short_name):
    """Reads in data from a target's .profile file.

    :param short_name: the short name of the target
    :type short_name: str
    :return: dictionary of the profile data
    :rtype: dict
    """

    pines_path = pat.utils.pines_dir_check()
    profile_path = pines_path/('Objects/'+short_name+'/'+short_name.replace(' ','').lower()+'.profile')

    if not os.path.exists(profile_path):
        print('{} does not exist, creating profile file with default params!'.format(profile_path.name))
        pines_path = pat.utils.pines_dir_check()
        red_path = pines_path/('Objects/'+str(profile_path).split('/')[-2]+'/reduced/')
        red_images = np.array(natsorted(glob(str(red_path)+'/*.fits')))
        source_detect_image = red_images[40].split('/')[-1]
        exclude_lower_left = 'False'
        dimness_tolerance = '0.33'
        brightness_tolerance = '3.0'
        distance_from_target = '900'
        non_linear_limit = '3750'
        edge_tolerance = '40'
        guess_position_x = '700'
        guess_position_y = '382'

        line = 'source_detect_image  = {}\nexclude_lower_left   = {}\ndimness_tolerance    = {}\nbrightness_tolerance = {}\ndistance_from_target = {}\nnon_linear_limit     = {}\nedge_tolerance       = {}\nguess_position_x     = {}\nguess_position_y     = {}'.format(source_detect_image, exclude_lower_left, dimness_tolerance, brightness_tolerance, distance_from_target, non_linear_limit, edge_tolerance, guess_position_x, guess_position_y)

        f = open(profile_path, 'x')
        with open(profile_path, 'w') as f:
            f.writelines(line)

    df = pd.read_csv(profile_path, sep='=', header=None).transpose()
    new_header = df.iloc[0]
    df = df[1:]
    df.columns = new_header
    df.columns = df.columns.str.strip()

    output_dict = {'source_detect_image':df['source_detect_image'].iloc[0].strip(),
                'exclude_lower_left':bool(distutils.util.strtobool(df['exclude_lower_left'].iloc[0].strip())),
                'dimness_tolerance':float(df['dimness_tolerance'].iloc[0]),
                'brightness_tolerance':float(df['brightness_tolerance'].iloc[0]),
                'distance_from_target':float(df['distance_from_target'].iloc[0]),
                'non_linear_limit':float(df['non_linear_limit'].iloc[0]),
                'edge_tolerance':float(df['edge_tolerance'].iloc[0]),
                'guess_position_x':float(df['guess_position_x'].iloc[0]),
                'guess_position_y':float(df['guess_position_y'].iloc[0])}
                
    return output_dict