import pdb 
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.data.get_master_synthetic_image import get_master_synthetic_image
from pines_analysis_toolkit.utils.pines_login import pines_login
from astropy.io import fits 
from pines_analysis_toolkit.utils.quick_plot import quick_plot as qp 
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from photutils.background import Background2D
from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader
from pines_analysis_toolkit.photometry.detect_sources import detect_sources
import numpy as np 
import matplotlib.pyplot as plt
from pines_analysis_toolkit.observing.synthetic_image_maker import synthetic_image_maker
from scipy import signal
from astropy.modeling import models, fitting
from natsort import natsorted 
from glob import glob 
import pandas as pd 
from pines_analysis_toolkit.observing.shift_measurer import shift_measurer

pd.options.mode.chained_assignment = None  # default='warn'
target = 'SIMP0136'
log_path = '/Users/tamburo/Documents/PINES_analysis_toolkit/Logs/20151110_log.txt'
log = pines_log_reader(log_path) #Get log shifts
files = np.array(natsorted(glob('/Users/tamburo/Documents/PINES_analysis_toolkit/Objects/SIMP0136/reduced/20151110*.fits'))) #Get files to measure new shifts with
pdb.set_trace()
#files = files[63:]
#files = np.array(natsorted(glob('/Users/tamburo/Documents/PINES_analysis_toolkit/Objects/2MASS 0014-0838/reduced/20201001.744_red.fits'))) #Get files to measure new shifts with
x_diffs = np.zeros(len(files))
y_diffs = np.zeros(len(files))
for i in range(len(files)):
    filename = files[i].split('/')[-1]
    
    log_ind = np.where(log['Filename'] == filename.split('_')[0]+'.fits')[0][0]
    log_x_shift = float(log['X shift'][log_ind])
    log_y_shift = float(log['Y shift'][log_ind])
    measured_x_shift, measured_y_shift = shift_measurer(target, filename)

    #Make sure the measured shifts are real values. 
    if np.isnan(measured_x_shift) or np.isnan(measured_y_shift):
        raise RuntimeError('Found nans for shifts!')
        pdb.set_trace()
        
    x_diff = abs(log_x_shift - measured_x_shift)
    y_diff = abs(log_y_shift - measured_y_shift)
    x_diffs[i] = x_diff
    y_diffs[i] = y_diff
    print('{}, {} of {}.'.format(filename, i+1, len(files)))
    print('x diff: {:3.2f}, y diff: {:3.2f}'.format(x_diff, y_diff))
    print('measured x shift: {:4.1f}, measured y shift: {:4.1f}'.format(measured_x_shift, measured_y_shift))
    print('')

    #Overwrite the telescope's logged shifts with the new measurements. 
    log['X shift'][log_ind] = str(np.round(measured_x_shift,1))
    log['Y shift'][log_ind] = str(np.round(measured_y_shift,1))


#Write out the new log. 
log.to_csv(log_path, index=0)
