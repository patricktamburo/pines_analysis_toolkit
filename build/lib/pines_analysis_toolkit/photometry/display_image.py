from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
import pdb

def display_image(file):
    data = fits.open(file)[0].data
    avg,med,std = sigma_clipped_stats(data,sigma=3)
    fig,ax = plt.subplots(figsize=(6,5))
    plt.ion()
    im = ax.imshow(data,origin='lower',vmin=med,vmax=med+4*std)
    plt.colorbar(im,label='Counts')
    plt.show()
    pdb.set_trace()