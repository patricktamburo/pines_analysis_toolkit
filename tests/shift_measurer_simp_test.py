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
pd.options.mode.chained_assignment = None  # Suppress some useless warnings. 
from pines_analysis_toolkit.observing.shift_measurer import shift_measurer
from pines_analysis_toolkit.utils import pines_dir_check 

#-----------------------------------------------------------------------------------------------------------
target = 'SIMP0136'
date = '20151110'
#-----------------------------------------------------------------------------------------------------------


pines_path = pines_dir_check()
short_name = short_name_creator(target)
log_path = pines_path/('Logs/'+date+'_log.txt')
reduced_path = pines_path/('Objects/'+short_name+'/reduced/')
files = np.array(natsorted(glob(str(reduced_path)+'/'+date+'*.fits'))) #Get files to measure new shifts with
log = pines_log_reader(log_path) #Get log shifts

for i in range(len(files)):
    filename = files[i].split('/')[-1]
    #Figure out which file you're looking at and its position in the log. 
    log_ind = np.where(log['Filename'] == filename.split('_')[0]+'.fits')[0][0]

    #Measure the shifts. 
    measured_x_shift, measured_y_shift = shift_measurer(target, filename, short_name)

    #Make sure the measured shifts are real values. 
    if np.isnan(measured_x_shift) or np.isnan(measured_y_shift):
        raise RuntimeError('Found nans for shifts!')
        pdb.set_trace()
        
    print('{}, {} of {}.'.format(filename, i+1, len(files)))
    print('measured x shift: {:4.1f}, measured y shift: {:4.1f}'.format(measured_x_shift, measured_y_shift))
    print('')

    #Overwrite the telescope's logged shifts with the new measurements. 
    log['X shift'][log_ind] = str(np.round(measured_x_shift,1))
    log['Y shift'][log_ind] = str(np.round(measured_y_shift,1))


#Write out the new log. 
log.to_csv(log_path, index=0)
