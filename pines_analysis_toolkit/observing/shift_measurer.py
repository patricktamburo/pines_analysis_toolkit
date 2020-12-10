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

def shift_measurer(target, image_name):
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
       
    

    pines_path = pines_dir_check()
    short_name = short_name_creator(target)

    synthetic_filename = target.replace(' ', '')+'_master_synthetic.fits'
    synthetic_path = pines_path/('Calibrations/Master Synthetic Images/'+synthetic_filename)

    image_path = pines_path/('Objects/'+short_name+'/reduced/'+image_name)
    #Check for the appropriate master synthetic image on disk. If it's not there, get from PINES server. 
    if not synthetic_path.exists():
        sftp = pines_login()
        get_master_synthetic_image(sftp, target)
        sftp.close()

    master_synthetic_image = fits.open(synthetic_path)[0].data

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
    seeing = log['X seeing'][ind]

    #Find sources in the image. 
    sources = detect_sources(image_path, seeing, edge_tolerance=12, thresh=4.0)
    
    qp(check_image)
    sort_inds = np.argsort(np.array(sources['aperture_sum']))[::-1]
    source_x = np.array(sources['xcenter'])[sort_inds[15:30]]
    source_y = np.array(sources['ycenter'])[sort_inds[15:30]]
    plt.plot(source_x, source_y, 'rx')
    pdb.set_trace()
    
    check_synthetic_image = synthetic_image_maker(source_x, source_y)

    #Measure the shift between the synthetic images.
    corr = signal.fftconvolve(master_synthetic_image,check_synthetic_image[::-1,::-1])

    (x_shift,y_shift) = corr_shift_determination(corr)
    print('(X shift, Y shift): ({:3.1f}, {:3.1f})'.format(x_shift, -y_shift))
    print('')

if __name__ == '__main__':
    filename = '20201205.346_red.fits'
    shift_measurer('2MASS J04574602-0207179', filename)
