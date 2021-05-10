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
from astropy.modeling import models, fitting
from astropy.stats import sigma_clipped_stats, sigma_clip

def log_updater(date, sftp, shift_tolerance=30, upload=False):
    '''
    Authors:
		Patrick Tamburo, Boston University, January 2021
	Purpose:
        Updates x_shift and y_shift measurements from a PINES log. These shifts are measured using *full* resolution images, while at the telescope,
        we use *half* resolution images (to save time between exposures). By measuring on full-res images, we get more accurate shifts, which allows 
        us to determine centroids more easily.
	Inputs:
        date (str): the UT date of the log whose shifts you want to update in YYYYMMDD format, e.g. '20151110'
        sftp (pysftp connection): sftp connection to the PINES server
        shift_tolerance (float): the maximum distance an x/y shift can be before shifts will be flagged as poor quality. 
        upload (bool): whether or not to push the updated log to the PINES server (only admins can do this)
    Outputs:
		Writes updated log file to disk. 
	TODO:
        Re-measure seeing?
    FIXME:
    '''
    def tie_sigma(model):
        return model.x_stddev_1

    def guide_star_seeing(subframe):
        # subframe = subframe - np.median(subframe)
        subframe = subframe - np.percentile(subframe,5)
        sub_frame_l = int(np.shape(subframe)[0])
        y, x = np.mgrid[:sub_frame_l, :sub_frame_l]

        # Fit with constant, bounds, tied x and y sigmas and outlier rejection:
        gaussian_init = models.Const2D(0.0) + models.Gaussian2D(subframe[int(sub_frame_l/2),int(sub_frame_l/2)],int(sub_frame_l/2),int(sub_frame_l/2),8/2.355,8/2.355,0)
        gaussian_init.x_stddev_1.min = 1.0/2.355
        gaussian_init.x_stddev_1.max = 20.0/2.355
        gaussian_init.y_stddev_1.min = 1.0/2.355
        gaussian_init.y_stddev_1.max = 20.0/2.355
        gaussian_init.y_stddev_1.tied = tie_sigma
        gaussian_init.theta_1.fixed = True
        fit_gauss = fitting.FittingWithOutlierRemoval(fitting.LevMarLSQFitter(),sigma_clip,niter=3,sigma=3.0)
        # gaussian, mask = fit_gauss(gaussian_init, x, y, subframe)
        gain = 8.21 #e per ADU
        read_noise = 2.43 #ADU
        weights = gain / np.sqrt(np.absolute(subframe)*gain + (read_noise*gain)**2) #1/sigma for each pixel
        gaussian, mask = fit_gauss(gaussian_init, x, y, subframe, weights)
        fwhm_x = 2.355*gaussian.x_stddev_1.value
        fwhm_y = 2.355*gaussian.y_stddev_1.value

        x_seeing = fwhm_x * 0.579
        y_seeing = fwhm_y * 0.579
        return(x_seeing,y_seeing)

    pines_path = pines_dir_check()
    log_path = pines_path/('Logs/'+date+'_log.txt')

    #Begin by checking filenames, making sure they're in sequential order, and that there is only one entry for each. 
    log_out_of_order_fixer(log_path, sftp)
    
    log = pines_log_reader(log_path) #Get telescope log shifts.
    myfile = open(log_path, 'r')
    lines = myfile.readlines()
    myfile.close()

    #Now loop over all files in the log, measure shifts in each file and update the line in the log. 
    for i in range(len(log)):
        if (log['Target'][i].lower() != 'flat') & (log['Target'][i].lower() != 'skyflat') & (log['Target'][i].lower() != 'supersky') & (log['Target'][i].lower() != 'dark') & (log['Target'][i].lower() != 'bias') & (log['Target'][i].lower() != 'dummy') & (log['Post-processing flag'][i] != 1):
            filename = log['Filename'][i].split('.fits')[0]+'_red.fits'
            target = log['Target'][i]
            short_name = short_name_creator(target)
            image_path = pines_path/('Objects/'+short_name+'/reduced/'+filename)

            #Figure out which file you're looking at and its position in the log. 
            log_ind = np.where(log['Filename'] == filename.split('_')[0]+'.fits')[0][0]

            #Measure the shifts and get positions of targets.
            (measured_x_shift, measured_y_shift, source_x, source_y, check_image) = shift_measurer(target, filename, sftp)

            if (abs(measured_x_shift) > shift_tolerance) or (abs(measured_y_shift) > shift_tolerance):
                print('Shift greater than {} pixels measured for {} in {}.'.format(shift_tolerance, short_name, image_path.name))
                print('Inspect manually.')
                shift_quality_flag = 1
            elif np.isnan(measured_x_shift) or np.isnan(measured_y_shift):
                raise RuntimeError('Found nans for shifts!')
                shift_quality_flag = 1
            else:
                shift_quality_flag = 0
            
            
            #Measure the seeing. 
            guide_star_cut = np.where((source_x > 50) & (source_x < 975) & (source_y > 50) & (source_y < 975))[0]
            if len(guide_star_cut) != 0:
                x_seeing_array = []
                y_seeing_array = []
                for guide_star_ind in guide_star_cut:
                    guide_star_x_int = int(source_x[guide_star_ind])
                    guide_star_y_int = int(source_y[guide_star_ind])
                    guide_star_subframe = check_image[guide_star_y_int-15:guide_star_y_int+15,guide_star_x_int-15:guide_star_x_int+15]
                    (x_seeing,y_seeing) = guide_star_seeing(guide_star_subframe)
                    #Cut unrealistic values/saturated stars. 
                    if x_seeing > 1.2 and x_seeing < 7.0:
                        x_seeing_array.append(x_seeing)
                        y_seeing_array.append(y_seeing)
                x_seeing = np.nanmedian(x_seeing_array)
                y_seeing = np.nanmedian(y_seeing_array)
            else:
                #Default to the average PINES value if no sources were found for guiding. 
                x_seeing = 2.6
                y_seeing = 2.6

            print('Log line {} of {}.'.format(i+1, len(log)))
            print('Measured x shift: {:4.1f}, measured y shift: {:4.1f}'.format(measured_x_shift, measured_y_shift))
            print('Measured seeing: {:4.1f}'.format(x_seeing))
            print('')

            #Overwrite the telescope's logged shifts and seeing values with the new measurements. 
            log['X shift'][log_ind] = str(np.round(measured_x_shift,1))
            log['Y shift'][log_ind] = str(np.round(measured_y_shift,1))
            log['X seeing'][log_ind] = str(np.round(x_seeing, 1))
            log['Y seeing'][log_ind] = str(np.round(y_seeing, 1))

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
            post_processing_flag = 1
            #Generate line of log text following the PINES telescope log format. 
            log_text = pines_logging(filename, log_date, target_name, filter_name, exptime, airmass, x_shift, y_shift, x_seeing, y_seeing, post_processing_flag, shift_quality_flag)

            #Overwrite the line with the new shifts.
            line_ind = log_ind + 1
            lines[line_ind] = log_text

            #Update the log on disk.
            with open(log_path, 'w') as f:
                for line in lines:
                    f.write(line)

        elif (log['Post-processing flag'][i] == 1):
            print('File already post-processed, skipping. {} of {}'.format(i+1, len(log)))
        else:
            print('File not a science target, skipping. {} of {}.'.format(i+1, len(log)))

    if upload:
        sftp.chdir('/data/logs/')
        print('Uploading to /data/logs/{}_log.txt.'.format(date))
        sftp.put(log_path,date+'_log.txt')
    return 