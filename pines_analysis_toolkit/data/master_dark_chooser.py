from datetime import datetime
import numpy as np
from astropy.io import fits


def master_dark_chooser(dark_path, header):
    """Identifies the closest dark in time to the image you want reduce, that has the correct exposure time.

    :param dark_path:  path to the ~.../PINES_analysis_toolkit/Calibrations/darks directory. 
    :type dark_path: pathlib.PosixPath
    :param header: header of the image you want to reduce
    :type header: astropy.io fits header
    :return: path to and name of the selected master dark 
    :rtype: pathlib.PosixPath, str
    """

    exptime = header['EXPTIME']
    obs_date = datetime.strptime(header['DATE-OBS'].split('T')[0].replace('-',''),'%Y%m%d')
    possible_darks = [x for x in (dark_path/'Master Darks').glob('*'+str(exptime)+'*.fits')]
    if (len(possible_darks) == 0): 
        print('ERROR: Could not find any suitable darks to reduce {}.'.format(header['FILENAME']))
        return
    else: 
        possible_dark_dates = [datetime.strptime(i.name.split('_')[-1].split('.')[0],'%Y%m%d') for i in possible_darks]
        dark_date_distances = [abs(possible_dark_dates[i]-obs_date) for i in range(len(possible_dark_dates))]
        dark_ind = np.where(np.array(dark_date_distances) == min(np.array(dark_date_distances)))[0][0]
        master_dark = fits.open(possible_darks[dark_ind])[0].data	
        master_dark_name = possible_darks[dark_ind].name
        return master_dark, master_dark_name