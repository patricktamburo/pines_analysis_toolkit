import pdb 
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.data.get_master_synthetic_image import get_master_synthetic_image
from pines_analysis_toolkit.utils.pines_login import pines_login
from pines_analysis_toolkit.data.reduce import reduce
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

def shift_measurer(target, image_name, sftp, num_sources=15, closeness_tolerance=10.):
    """Measure shifts between an image and the master synthetic image for a target. 

    :param target: target's full 2MASS name 
    :type target: str
    :param image_name: name of the image whose shifts you want to measure
    :type image_name: str
    :param sftp: sftp connection to the PINES server
    :type sftp: pysftp.connection
    :param num_sources: number of sources to use for measuring shift, defaults to 15
    :type num_sources: int, optional
    :param closeness_tolerance: closest targets can be in pixels and still be considered for correlation sources, defaults to 10.
    :type closeness_tolerance: float, optional
    :return: x/y shift, source x/y, check_image
    :rtype: [type]
    """
    '''Authors:
            Patrick Tamburo, Boston University, June 2020
        Purpose:
            Measures shifts between a particular image (a "check" image) and the master image for the field.
        Inputs:
            target (str): The target's full 2MASS name
            image_name (str): The name of the reduced image whose shift you want to measure. 
            num_sources (int): The maximum number of sources you would like to use to measure shifts. 
            closeness_tolerance (float): The closest two identified sources can be in pixels before one is thrown out. 
        Outputs:

        TODO:
            Grab logs automatically? 
            Flag bad centroids?
    '''

    def corr_shift_determination(corr):
        #Measure shift between the check and master images by fitting a 2D gaussian to corr. This gives sub-pixel accuracy. 
        y_max, x_max = np.unravel_index(np.argmax(corr), corr.shape) #Find the pixel with highest correlation, then use this as estimate for gaussian fit.
        y, x = np.mgrid[y_max-10:y_max+10, x_max-10:x_max+10]
        corr_cut =  corr[y,x]
        gaussian_init = models.Gaussian2D(np.max(corr_cut),x_max,y_max,8/2.355,8/2.355,0)
        fit_gauss = fitting.LevMarLSQFitter()
        gaussian = fit_gauss(gaussian_init, x, y,corr_cut)
        fit_x = gaussian.x_mean.value
        fit_y = gaussian.y_mean.value

        x_shift = fit_x - 1024
        y_shift = fit_y - 1024 
        return(x_shift,y_shift)
       
    #image_name = '20201005.354_red.fits'

    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    synthetic_filename = target.replace(' ', '')+'_master_synthetic.fits'
    synthetic_path = pines_path/('Calibrations/Master Synthetic Images/'+synthetic_filename)

    image_path = pines_path/('Objects/'+short_name+'/reduced/'+image_name)
    #Check for the appropriate master synthetic image on disk. If it's not there, get from PINES server. 
    if not synthetic_path.exists():
        get_master_synthetic_image(sftp, target)

    master_synthetic_image = fits.open(synthetic_path)[0].data
    
    #If the reduced image doesn't exist, download raw from PINES server and reduce it. 
    if not image_path.exists():
        missing_filename = image_path.name.split('_')[0]+'.fits'
        sftp.chdir('/data/raw/mimir/')
        runs = sftp.listdir()
        run_guess = missing_filename[0:6]
        ind = np.where(np.array(runs) == run_guess)[0][0]

        if ind + 1 == len(runs):
            inds = np.arange(ind-1, ind+1)
            runs = np.array(runs)[inds]
            runs = [runs[1], runs[0]]
        else:
            inds = np.arange(ind-1, ind+2)
            runs = np.array(runs)[inds]
            runs = [runs[1], runs[0], runs[2]] 

        for jj in range(len(runs)):
            nights = sftp.listdir('/data/raw/mimir/'+runs[jj])
            nights = [i for i in nights if i[0] == '2']
            for kk in range(len(nights)):
                server_files = sftp.listdir('/data/raw/mimir/'+runs[jj]+'/'+nights[kk])
                if missing_filename in server_files:
                    print('File found: {}'.format('/data/raw/mimir/'+runs[jj]+'/'+nights[kk]+'/'+missing_filename))
                    found = True
                    #Download the file and grab relevant parameters from the header. 
                    sftp.get('/data/raw/mimir/'+runs[jj]+'/'+nights[kk]+'/'+missing_filename, pines_path/('Objects/'+short_name+'/raw/'+missing_filename))
                    break
        reduce(target)

    #Read in the check image and interpolate/background subtract to make source detection easier. 
    check_image = fits.open(image_path)[0].data
    check_image = interpolate_replace_nans(check_image, kernel=Gaussian2DKernel(x_stddev=0.5))
    bg_2d = Background2D(check_image, box_size=64)
    check_image = check_image - bg_2d.background

    #Read in the log and figure out what seeing FWHM to use. 
    date = image_name.split('.')[0]
    log = pines_log_reader(pines_path/('Logs/'+date+'_log.txt'))
    raw_filename = image_name.split('_')[0]+'.fits'
    ind = np.where(log['Filename'] == raw_filename)[0][0]
    seeing = float(log['X seeing'][ind])

    if (seeing <= 1.0) or (seeing >= 7.0) or (np.isnan(seeing)):
        if ind >= 5:
            seeing = np.nanmedian(np.array(log['X seeing'][ind-5:ind], dtype='float'))
        else:
            seeing = 2.6
    


    #Find sources in the image. 
    sources = detect_sources(image_path, seeing, edge_tolerance=10, thresh=3.5)
    
    #Comb the returned sources and cut any that are too close to another one. 
    #This can happen if the actual seeing differs from that recorded in the log. 
    bad_inds = []
    for i in range(len(sources)):
        if i not in bad_inds:
            x = sources['xcenter'][i]
            y = sources['ycenter'][i]
            dists = np.array(np.sqrt((sources['xcenter']-x)**2 + (sources['ycenter']-y)**2))
            duplicate_inds = np.where((dists < closeness_tolerance) & (dists != 0))[0]
            if len(duplicate_inds) >= 1:
                bad_inds.extend(duplicate_inds)
    

    ap_sum = np.array(sources['aperture_sum'])
    source_x = np.array(sources['xcenter'])
    source_y = np.array(sources['ycenter'])

    ap_sum = np.delete(ap_sum, bad_inds)
    source_x = np.delete(source_x, bad_inds)
    source_y = np.delete(source_y, bad_inds)

    sort_inds = np.argsort(ap_sum)[::-1]
    source_x = source_x[sort_inds[0:num_sources]]
    source_y = source_y[sort_inds[0:num_sources]]

    # if raw_filename == '20201005.354.fits':
    #     qp(check_image)
    #     plt.plot(source_x, source_y, 'rx')
    #     pdb.set_trace()
    
    check_synthetic_image = synthetic_image_maker(source_x, source_y)

    #Measure the shift between the synthetic images.
    corr = signal.fftconvolve(master_synthetic_image,check_synthetic_image[::-1,::-1])

    (x_shift,y_shift) = corr_shift_determination(corr)
    #print('(X shift, Y shift): ({:3.1f}, {:3.1f})'.format(x_shift, -y_shift))
    #print('')



    return (x_shift, -y_shift, source_x, source_y, check_image)