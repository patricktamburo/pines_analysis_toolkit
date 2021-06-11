import pdb
import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ImageNormalize, ZScaleInterval
from astropy.convolution import interpolate_replace_nans, Gaussian2DKernel

def quick_plot(array, interp=False):

    if interp:
        array = interpolate_replace_nans(array, kernel=Gaussian2DKernel(x_stddev=0.25))

    plt.ion()
    plt.figure()
    norm = ImageNormalize(array, interval=ZScaleInterval())
    im = plt.imshow(array, origin='lower', norm=norm)
    cb = plt.colorbar(im)
    plt.tight_layout()