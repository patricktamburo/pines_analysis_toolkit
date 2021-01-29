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
from pines_analysis_toolkit.observing.pines_logging import pines_logging
from pines_analysis_toolkit.utils.pines_login import pines_login
from pines_analysis_toolkit.observing.log_out_of_order_fixer import log_out_of_order_fixer

def log_updater(target, date, upload=False):
    '''
    Authors:
		Patrick Tamburo, Boston University, January 2021
	Purpose:
        Updates x_shift and y_shift measurements from a PINES log. These shifts are measured using *full* resolution images, while at the telescope,
        we use *half* resolution images (to save time between exposures). By measuring on full-res images, we get more accurate shifts, which allows 
        us to determine centroids more easily.
	Inputs:
        target (str): a targets 'long' 2MASS name, e.g. '2MASS J01234567+012345678'
        date (str): the UT date of the log whose shifts you want to update in YYYYMMDD format, e.g. '20151110'
        upload (bool): whether or not to push the updated log to the PINES server (only admins can do this)
    Outputs:
		Writes updated log file to disk. 
	TODO:
        Re-measure seeing?
    FIXME:
    '''

    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    log_path = pines_path/('Logs/'+date+'_log.txt')
    reduced_path = pines_path/('Objects/'+short_name+'/reduced/')
    files = np.array(natsorted(glob(str(reduced_path)+'/'+date+'*.fits'))) #Get files to measure new shifts with

    #Begin by checking filenames, making sure they're in sequential order, and that there is only one entry for each. 
    log_out_of_order_fixer(log_path)
    
    log = pines_log_reader(log_path) #Get telescope log shifts.
    myfile = open(log_path, 'r')
    lines = myfile.readlines()
    myfile.close()

    #Now loop over all files for this target in the log, measure shifts in each file and update the line in the log. 
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

        #Grab entries for log line.
        filename = log['Filename'][log_ind]
        log_date = log['Date'][log_ind]
        target_name = log['Target'][log_ind]
        filter_name = log['Filt.'][log_ind]
        exptime = log['Exptime'][log_ind]
        airmass = log['Airmass'][log_ind]
        x_shift = log['X shift'][log_ind]
        y_shift = log['Y shift'][log_ind]
        x_seeing = log['X seeing'][log_ind]
        y_seeing = log['Y seeing'][log_ind]
        
        #Generate line of log text following the PINES telescope log format. 
        log_text = pines_logging(filename, log_date, target_name, filter_name, exptime, airmass, x_shift, y_shift, x_seeing, y_seeing)

        #Overwrite the line with the new shifts.
        line_ind = log_ind + 1
        lines[line_ind] = log_text

        #Update the log on disk.
        with open(log_path, 'w') as f:
            for line in lines:
                f.write(line)

    if upload:
        sftp = pines_login()
        sftp.chdir('/data/logs/')
        print('Uploading to /data/logs/{}_log.txt.'.format(date))
        sftp.put(log_path,date+'_log.txt')
        sftp.close()
    return 