import pdb
import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ImageNormalize, ZScaleInterval
def quick_plot(array):
    plt.ion()
    plt.figure()
    norm = ImageNormalize(array, interval=ZScaleInterval())
    im = plt.imshow(array, origin='lower', norm=norm)
    cb = plt.colorbar(im)
    plt.tight_layout()