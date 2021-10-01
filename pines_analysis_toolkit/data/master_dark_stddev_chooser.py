from datetime import datetime
import numpy as np 
from astropy.io import fits

def master_dark_stddev_chooser(dark_std_path, header):
    """Identifies the closest master_dark_stddev in time to the image you want reduce, that has the correct exposure time.

    :param dark_std_path:  path to the ~.../PINES_analysis_toolkit/Calibrations/dark_std directory. 
    :type dark_std_path: pathlib.PosixPath
    :param header: header of the image you want to reduce
    :type header: astropy.io fits header
    :return: path to and name of the selected master dark stddev file
    :rtype: pathlib.PosixPath, str
    """

    exptime = header['EXPTIME']
    obs_date = datetime.strptime(header['DATE-OBS'].split('T')[0].replace('-',''),'%Y%m%d')
    possible_dark_stds = [x for x in dark_std_path.glob('*'+str(exptime)+'*.fits')]
    if (len(possible_dark_stds)) == 0:
        print('ERROR: Could not find any suitable master_dark_stddev images to get read noise/dark current measurement for {}.'.format(header['FILENAME']))
        return
    else:
        possible_dark_std_dates = [datetime.strptime(i.name.split('_')[-1].split('.')[0],'%Y%m%d') for i in possible_dark_stds]
        dark_std_date_distances = [abs(possible_dark_std_dates[i]-obs_date) for i in range(len(possible_dark_std_dates))]
        dark_std_ind = np.where(np.array(dark_std_date_distances) == min(np.array(dark_std_date_distances)))[0][0]
        master_dark_std = fits.open(possible_dark_stds[dark_std_ind])[0].data
        master_dark_std_name = possible_dark_stds[dark_std_ind].name
        return master_dark_std