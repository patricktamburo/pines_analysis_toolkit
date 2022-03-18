import pines_analysis_toolkit as pat

from astropy.time import Time
from astropy import coordinates as coord
from astropy import units as u
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ImageNormalize, ZScaleInterval
from astropy.convolution import interpolate_replace_nans, Gaussian2DKernel

import pandas as pd
import pysftp
import getpass
import distutils
import os
from glob import glob
from natsort import natsorted
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from scipy.stats import sigmaclip
from datetime import datetime
import imageio

def bjd_tdb_to_jd_utc(bjd_tdb, ra, dec, location='lowell'):
    """Converts a Barycentric Julian Date in TDB timescale to Julian Date in UTC timescale.

    :param bjd_tdb: the barycentric julian date in tdb
    :type bjd_tdb: float
    :param ra: the ra of the target (e.g., '23:06:30.0')
    :type ra: str
    :param dec: the dec of the target (e.g, '-05:01:57')
    :type dec: str
    :param location: the site of observations, must match a location in the astropy .of_site json file. , defaults to 'lowell'
    :type location: str, optional
    :return: time in jd utc
    :rtype: float
    """
    if location == 'lowell':
        #Specify the location of lowell without needing to download the sites.json file.
        #The values used are the geocentric coordinates of lowell in meters, and I got them from:
        #   coord.EarthLocation.of_site('lowell')
        site = coord.EarthLocation.from_geocentric(-1918329.73705223*u.m, -4861253.21396165*u.m, 3647910.33904671*u.m)
    else:
        site = coord.EarthLocation.of_site(location)

    input_bjd_tdb = Time(bjd_tdb, format='jd', scale='tdb', location=site)
    target = coord.SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
    ltt_bary = input_bjd_tdb.light_travel_time(target)
    return (input_bjd_tdb.utc - ltt_bary).value


def get_source_names(df):
    """Grabs sources names from a PINES centroided sources dataframe.

    :param df: dataframe from which you want to grab source names
    :type df: pandas DataFrame
    :return: list of source names
    :rtype: list of str
    """
    df_keys = df.keys().str.strip().str.replace(' Corrected Flux Error', '').str.replace(' Image X', '').str.replace(' Image Y', '').str.replace(' Cutout X', '').str.replace(' Cutout Y', '').str.replace(' Centroid Warning', '').str.replace(' Corrected Flux', '').str.replace(
        ' ALC Weight', '').str.replace(' Background', '').str.replace(' Raw Flux', '').str.replace(' Bad Pixel Flag', '').str.replace(' Interpolation Flag', '').str.replace(' Flux', '').str.replace(' Flux Error', '').str.replace(' Error', '').str.replace(' Normalized', '')
    source_names = []
    for i in range(len(df_keys)):
        df_key = df_keys[i]
        if (df_key != 'Filename') and (df_key != 'Night Number') and (df_key != 'Block Number') and (df_key != 'Seeing') and (df_key != 'Airmass') and (df_key != 'Time BJD TDB') and (df_key != 'Time (JD UTC)') and (df_key != 'Time UT') and (df_key != 'Time JD UTC') and (df_key != 'Sigma Clip Flag') and (df_key != 'Block Number') and (df_key != 'ALC') and (df_key not in source_names):
            source_names.append(df_key)
    return source_names


def gif_maker(path='.', ext='.jpg', fps=30, title='animation'):
    """Creates a gif from a directory of images. 

    :param path: path to file directory containing images, defaults to '.'
    :type path: str, optional
    :param ext: extension of images, defaults to '.jpg'
    :type ext: str, optional
    :param fps: frames per second in the gif, defaults to 30
    :type fps: int, optional
    :param title: name of the output file +.gif, defaults to 'animation'
    :type title: str, optional
    """
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
    print('gif size: {} Mb'.format(np.round(
        os.stat(file_path/('animation/'+title+'.gif')).st_size / (1000**2), 1)))


def jd_utc_to_bjd_tdb(jd_utc, ra, dec, location='lowell'):
    """Converts Julian Date in UTC timescale to Barycentric Julian Date in TDB timescale. 

    :param jd_utc: julian date in utc
    :type jd_utc: float
    :param ra: the ra of the target as an hour angle (e.g., '23:06:30.0') or decimal degrees (e.g., 346.622013)
    :type ra: str or float
    :param dec: the dec of the target as a string (e.g, '-05:01:57') or decimal degrees (e.g., -5.041274)
    :type dec: str or float
    :param location: the site of observations, must match a location in the astropy .of_site json file, defaults to 'lowell'
    :type location: str, optional
    :return: time in bjd tdb
    :rtype: float
    """

    if type(ra) == str:
        ra_unit = u.hourangle
    elif type(ra) == np.float64:
        ra_unit = u.deg

    if location == 'lowell':
        #Specify the location of lowell without needing to download the sites.json file.
        #The values used are the geocentric coordinates of lowell in meters, and I got them from:
        #   coord.EarthLocation.of_site('lowell')
        site = coord.EarthLocation.from_geocentric(-1918329.73705223*u.m, -4861253.21396165*u.m, 3647910.33904671*u.m)
    else:
        site = coord.EarthLocation.of_site(location)

    input_jd_utc = Time(jd_utc, format='jd', scale='utc', location=site)
    target = coord.SkyCoord(ra, dec, unit=(ra_unit, u.deg), frame='icrs')
    ltt_bary = input_jd_utc.light_travel_time(target)
    return (input_jd_utc.tdb + ltt_bary).value


def object_directory_creator(pines_path, short_name):
    """Sets up the directory structure for an object.

    :param pines_path: path to the top-level PINES directory, e.g. ~/Documents/PINES_analysis_toolkit/
    :type pines_path: pathlib.PosixPath
    :param short_name: short name for the target
    :type short_name: str
    """
    os.mkdir(pines_path/('Objects/'+short_name))
    os.mkdir(pines_path/('Objects/'+short_name+'/analysis'))
    os.mkdir(pines_path/('Objects/'+short_name+'/analysis/diagnostic_plots'))
    os.mkdir(pines_path/('Objects/'+short_name+'/analysis/aper_phot_analysis'))
    os.mkdir(pines_path/('Objects/'+short_name+'/aper_phot'))
    os.mkdir(pines_path/('Objects/'+short_name+'/output'))
    os.mkdir(pines_path/('Objects/'+short_name+'/psf_phot'))
    os.mkdir(pines_path/('Objects/'+short_name+'/pwv'))
    os.mkdir(pines_path/('Objects/'+short_name+'/sources'))
    os.mkdir(pines_path/('Objects/'+short_name+'/raw'))
    os.mkdir(pines_path/('Objects/'+short_name+'/reduced'))
    return


def pines_dir_check():
    """Finds the local PINES_analysis_toolkit directory

    :return: path to PINES directory
    :rtype: pathlib.PosixPath
    """
    home_dir = Path(os.path.expanduser('~/Documents/'))
    default_pines_dir = Path(os.path.join(home_dir, 'PINES_analysis_toolkit/'))
    if os.path.exists(default_pines_dir):
        return default_pines_dir
    else:
        print('ERROR...I have not set this up to work for directories other than ~/Documents/PINES_analysis_toolkit.')
        # TODO: Config file for non-Documents path?
        return


def pines_log_reader(path):
    """Reads a PINES log or csv file into a dataframe

    :param path: path to the log.txt file or .csv file.
    :type path: pathlib.PosixPath
    :return: datframe of the log/csv file data
    :rtype: pandas DataFrame
    """
    df = pd.read_csv(path)

    # Remove trailing/leading spaces from column names
    df.columns = df.columns.str.lstrip()
    df.columns = df.columns.str.rstrip()

    # Remove our header comment idicator in the first column if it's there.
    if '#' in df.columns[0]:
        df.rename(
            columns={df.columns[0]: df.columns[0].replace('#', '')}, inplace=True)

    # Remove trailing and leading spaces from log entries.
    for key in df.keys():
        try:
            df[key] = df[key].str.strip()
        except:
            continue

    return df

def mimir_date_reader(header):
    date_obs = header['DATE-OBS']
    # Catch a case that can cause datetime strptime to crash; Mimir headers sometimes have DATE-OBS with seconds specified as 010.xx seconds, when it should be 10.xx seconds.
    if len(date_obs.split(':')[-1].split('.')[0]) == 3:
        date_obs = date_obs.split(
            ':')[0] + ':' + date_obs.split(':')[1] + ':' + date_obs.split(':')[-1][1:]

    if date_obs.split(':')[-1] == '60.00':
        date_obs = date_obs.split(
            ':')[0]+':'+str(int(date_obs.split(':')[1])+1)+':00.00'
    if '.' not in date_obs:
        time_fmt = '%Y-%m-%dT%H:%M:%S'
    else:
        time_fmt = '%Y-%m-%dT%H:%M:%S.%f'

    # Keep a try/except clause here in case other unknown DATE-OBS formats pop up.
    try:
        date = datetime.strptime(date_obs, time_fmt)
        return date
    except:
        print(
            'Header DATE-OBS format does not match the format code in strptime! Inspect/correct the DATE-OBS value.')
        breakpoint()


def pines_login():
    """Establishes an sftp connection with the PINES server.

    :return: sftp connection
    :rtype: pysftp connection
    """
    i = 3
    while i != 0:
        print('')
        username = input('Enter pines.bu.edu username: ')
        password = getpass.getpass(
            'Enter password for {}@pines.bu.edu: '.format(username))
        print('')
        try:
            sftp = pysftp.Connection(
                'pines.bu.edu', username=username, password=password)
            print('Login successful!\n')
            return sftp
        except:
            i -= 1
            if i == 1:
                verb = 'try'  # lol
            else:
                verb = 'tries'
            if i != 0:
                print('Login failed! {} {} remaining.'.format(i, verb))
            else:
                print('ERROR: login failed after {} tries.'.format(i))


def profile_reader(short_name, force_output_path=''):
    """Reads in data from a target's .profile file.

    :param short_name: the short name of the target
    :type short_name: str
    :param force_output_path: user-chosen directory to use in place of the default ~/Documents/PINES_analysis_toolkit directory for analysis, defaults to ''
    :type force_output_path: pathlib.PosixPath
    :return: dictionary of the profile data
    :rtype: dict
    """
    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()    

    profile_path = pines_path/('Objects/'+short_name+'/'+short_name.replace(' ', '').lower()+'.profile')


    if not os.path.exists(profile_path):
        print('{} does not exist, creating profile file with default params!'.format(
            profile_path.name))
        pines_path = pat.utils.pines_dir_check()
        red_path = pines_path/('Objects/'+str(profile_path).split('/')[-2]+'/reduced/')
        red_images = np.array(natsorted(glob(str(red_path)+'/*.fits')))
        source_detect_image = red_images[40].split('/')[-1]
        exclude_lower_left = 'False'
        dimness_tolerance = '0.40'
        brightness_tolerance = '3.0'
        distance_from_target = '900'
        non_linear_limit = '3750'
        edge_tolerance = '40'
        guess_position_x = '700'
        guess_position_y = '382'

        line = 'source_detect_image  = {}\nexclude_lower_left   = {}\ndimness_tolerance    = {}\nbrightness_tolerance = {}\ndistance_from_target = {}\nnon_linear_limit     = {}\nedge_tolerance       = {}\nguess_position_x     = {}\nguess_position_y     = {}'.format(
            source_detect_image, exclude_lower_left, dimness_tolerance, brightness_tolerance, distance_from_target, non_linear_limit, edge_tolerance, guess_position_x, guess_position_y)

        f = open(profile_path, 'x')
        with open(profile_path, 'w') as f:
            f.writelines(line)

    df = pd.read_csv(profile_path, sep='=', header=None).transpose()
    new_header = df.iloc[0]
    df = df[1:]
    df.columns = new_header
    df.columns = df.columns.str.strip()

    output_dict = {'source_detect_image': df['source_detect_image'].iloc[0].strip(),
                   'exclude_lower_left': bool(distutils.util.strtobool(df['exclude_lower_left'].iloc[0].strip())),
                   'dimness_tolerance': float(df['dimness_tolerance'].iloc[0]),
                   'brightness_tolerance': float(df['brightness_tolerance'].iloc[0]),
                   'distance_from_target': float(df['distance_from_target'].iloc[0]),
                   'non_linear_limit': float(df['non_linear_limit'].iloc[0]),
                   'edge_tolerance': float(df['edge_tolerance'].iloc[0]),
                   'guess_position_x': float(df['guess_position_x'].iloc[0]),
                   'guess_position_y': float(df['guess_position_y'].iloc[0])}

    return output_dict


def quick_plot(array, interp=False):
    """Creates a quick plot of a 2D array using a ZScaleInterval

    :param array: 2D image array
    :type array: numpy array
    :param interp: whether or not to interpolate any NaNs in the image, defaults to False
    :type interp: bool, optional
    """
    if interp:
        array = interpolate_replace_nans(
            array, kernel=Gaussian2DKernel(x_stddev=0.25))

    plt.ion()
    plt.figure()
    norm = ImageNormalize(array, interval=ZScaleInterval())
    im = plt.imshow(array, origin='lower', norm=norm)
    cb = plt.colorbar(im)
    plt.tight_layout()


def ref_counter():
    """Counts the number of reference stars that each object in your ~/PINES_analysis_toolkit/Objects/ directory. 
    Prints the mean of the distribution and +/- 1-sigma range. 
    """
    pines_path = pat.utils.pines_dir_check()
    objects_path = pines_path/('Objects/')
    targets = natsorted(os.listdir(objects_path))
    targets = np.array(targets[0:-1])
#    # Remove staring sdM target.
#    targets = np.delete(targets, np.where(targets == '2MASS 0438+4349')[0][0])
#    # Remove target that is saturated.
#    targets = np.delete(targets, np.where(targets == '2MASS 1333-0215')[0][0])
    num_refs = np.zeros(len(targets))
    for i in range(len(targets)):
        source_path = objects_path / \
            (targets[i]+'/sources/target_and_references_source_detection.csv')
        try:
            df = pd.read_csv(source_path)
            num_refs[i] = len(df) - 1  # Subtract the target.
        except:
            num_refs[i] = np.nan

    num_refs = num_refs[~np.isnan(num_refs)]
    v, l, h = sigmaclip(num_refs, 3, 3)
    res = np.percentile(v, [15.9, 50, 84.1])
    print('{}^+{}_-{}'.format(res[1], res[1]-res[0], res[2]-res[1]))

def short_name_creator(long_name):
    """Tries to create a 'short name' out of a full 2MASS identifier. 
    I.e., 2MASS J00001354+2554180 would be converted to 2MASS 0000+2554. 

    :param long_name: full 2MASS name of the target
    :type long_name: str
    :return: short name for the target 
    :rtype: str
    """
    name = long_name.replace(' ', '')

    if 'UT' in name:
        long_name = long_name.split('UT')[0].rstrip()
        name = name.split('UT')[0]

    if '2MASS' in name:
        # If the user already passed a 2MASS short name (i.e., 2MASS XXXX+YYYY), just return it.
        if len(name) == 14 or len(name) == 16 or len(name) == 17: #TODO: Fix this, stupid. 
            short_name = long_name

        # Otherwise, create the short name.
        else:
            if '+' in name:
                short_name = '2MASS ' + \
                    name.split('J')[1].split('+')[0][0:4]+'+' + \
                    name.split('J')[1].split('+')[1][0:4]
            elif '-' in name:
                short_name = '2MASS ' + \
                    name.split('J')[1].split('-')[0][0:4]+'-' + \
                    name.split('J')[1].split('-')[1][0:4]
    # If a name other than a 2MASS identifier was passed, just return the long name.
    else:
        print('')
        print('WARNING: long_name does not match 2MASS Jxxxxxxxx+xxxxxxx format, returning long_name.')
        short_name = long_name

    return short_name


def update_header(file_path, header_key, new_value):
    """Updates a fits header key with a chosen value

    :param file_path: path to the fits file
    :type file_path: pathlib.PosixPath
    :param header_key: header key you want to update (e.g., 'OBJECT', 'FILTNME2', etc.)
    :type header_key: str
    :param new_value: the value you want to use to update the header key 
    :type new_value: same type as header[header_key]
    """

    fits.setval(file_path, header_key, value=new_value)
