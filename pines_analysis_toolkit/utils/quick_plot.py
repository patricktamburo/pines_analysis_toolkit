import pdb
import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import sigma_clipped_stats

def quick_plot(array):
    plt.ion()
    plt.figure()
    avg, med, std = sigma_clipped_stats(array)
    im = plt.imshow(array, origin='lower', vmin=med, vmax=med+3*std)
    cb = plt.colorbar(im)
    plt.tight_layout()