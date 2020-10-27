import numpy as np
from astropy.stats import SigmaClip, sigma_clipped_stats
from photutils import Background2D, MedianBackground
import pdb 
import matplotlib.pyplot as plt

'''Authors:
		Patrick Tamburo, Boston University, August 2020
	Purpose:
		Removes large scale changes in the background of an image. See: https://photutils.readthedocs.io/en/stable/background.html
	Inputs:
        data (2d numpy array): Array of pixel values. 
        box_size (int, optional): Size of the boxes used to do the background estimation. 
	Outputs:
        data - bkg.background (2d numpy array): The image with background trends corrected. 
    TODO: 
'''

def bg_2d(data, box_size=50):
    sigma_clip = SigmaClip(sigma=3.)
    bkg_estimator = MedianBackground()
    bkg = Background2D(data, (box_size, box_size), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    # plt.ion()
    # fig, ax = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True, figsize=(16,7))
    # avg, med, std = sigma_clipped_stats(data)
    # ax[0].imshow(data, origin='lower', vmin=med, vmax=med+3*std)
    # ax[1].imshow(bkg.background, origin='lower', vmin=med, vmax=med+3*std)
    # avg, med, std = sigma_clipped_stats(data-bkg.background)
    # ax[2].imshow(data-bkg.background, origin='lower', vmin=med, vmax=med+3*std)
    # plt.tight_layout()
    # pdb.set_trace() 

    return data - bkg.background
