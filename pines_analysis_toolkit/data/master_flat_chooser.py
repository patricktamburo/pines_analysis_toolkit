from datetime import datetime
import numpy as np 
from astropy.io import fits

def master_flat_chooser(flats_path, header):
    """Identifies the closest master flat in time to the image you want reduce, that has the correct exposure time.

    :param flats_path:  path to the ~.../PINES_analysis_toolkit/Calibrations/Flats/Domeflats/ directory. 
    :type dark_path: pathlib.PosixPath
    :param header: header of the image you want to reduce
    :type header: astropy.io fits header
    :return: path to and name of the selected master flat
    :rtype: pathlib.PosixPath, str
    """
    band = header['FILTNME2']
    obs_date = datetime.strptime(header['DATE-OBS'].split('T')[0].replace('-',''),'%Y%m%d')
    possible_flats = [x for x in (flats_path/(band+'/Master Flats')).glob('*.fits')]
    if (len(possible_flats)) == 0:
        print('ERROR: Could not find any suitable flats to reduce {}.'.format(header['FILENAME']))
        return
    else:
        possible_flat_dates = [datetime.strptime(i.name.split('_')[-1].split('.')[0],'%Y%m%d') for i in possible_flats]
        flat_date_distances = [abs(possible_flat_dates[i]-obs_date) for i in range(len(possible_flat_dates))]
        flat_ind = np.where(np.array(flat_date_distances) == min(np.array(flat_date_distances)))[0][0]
        master_flat = fits.open(possible_flats[flat_ind])[0].data
        master_flat_name = possible_flats[flat_ind].name
        return master_flat, master_flat_name