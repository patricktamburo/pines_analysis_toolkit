import pdb
import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import sigma_clipped_stats

def quick_plot(array):
    plt.figure()
    avg, med, std = sigma_clipped_stats(array)
    plt.imshow(array, origin='lower', vmin=med, vmax=med+3*std)