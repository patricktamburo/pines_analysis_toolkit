from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
import matplotlib.pyplot as plt
import numpy as np
import pdb
import pandas as pd
from astropy.io import fits
import natsort
import datetime
import math
import scipy.ndimage

'''Authors:
		Patrick Tamburo, Boston University, September 2020
	Purpose:
        Measures drift rate in x, y and sqrt(x^2+y^2) for the target given a time series of centroids.
	Inputs:
        target (str): The target's full 2MASS name
        sources (pandas dataframe): List of source names, x and y positions. 
        plots (bool, optional): Whether or not to output plots showing centroid positions. Images output to sources directory within the object directory. 
        restore (bool, optional): Whether or not to restore centroider output that already exists. 
    Outputs:
		centroid_df (pandas DataFrame): X and Y centroid positions for each source. 
	TODO:
        Grab logs automatically? 
        Flag bad centroids?
'''

def drift_rate(target):
    plt.ion()

    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    source_centroids = pd.read_csv(pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv'))
    source_centroids.columns = source_centroids.columns.str.strip()

    x = source_centroids[short_name+' X']
    y = source_centroids[short_name+' Y']
    z = np.sqrt(x**2 + y**2)

    reduced_images = natsort.natsorted(list((pines_path/('Objects/'+short_name+'/reduced/')).rglob('*.fits')))

    if len(reduced_images) != len(x):
        raise RuntimeError('Number of reduced images does not match number of centroids.')
    
    dates = []
    for i in range(len(reduced_images)):
        file = reduced_images[i]
        date_obs = fits.open(file)[0].header['DATE-OBS']
        dates.append(datetime.datetime.strptime(date_obs, '%Y-%m-%dT%H:%M:%S.%f'))
    dates = np.array(dates)
    unique_days = natsort.natsorted(list({(i.day,i.month,i.year) for i in dates}))
    days = [(dates[i].day, dates[i].month, dates[i].year) for i in range(len(dates))]

    colors = ['tab:blue', 'tab:orange', 'tab:green']
    fig, ax = plt.subplots(nrows=3, ncols=len(unique_days), figsize=(17,8), sharey=True)
    for i in range(len(unique_days)):
        day = unique_days[i]
        truth_arr = [days[i] == day for i in range(len(days))]
        date_locs = np.where(truth_arr)[0]
        ax[i,0].plot(dates[date_locs], np.gradient(x[date_locs]), marker='o', ls='', color=colors[0], alpha=0.3)
        smoothed = scipy.ndimage.filters.uniform_filter(np.gradient(x[date_locs]),size=5)
        ax[i,0].plot(dates[date_locs], smoothed, lw=2, color=colors[0], zorder=0)
        for label in ax[i,0].get_xticklabels():
            label.set_rotation(30)
            label.set_horizontalalignment('right')

        ax[i,1].plot(dates[date_locs], np.gradient(y[date_locs]), marker='o', ls='', color=colors[1], alpha=0.3)
        smoothed = scipy.ndimage.filters.uniform_filter(np.gradient(y[date_locs]),size=5)
        ax[i,1].plot(dates[date_locs], smoothed,   lw=2, color=colors[1], zorder=0)
        for label in ax[i,1].get_xticklabels():
            label.set_rotation(30)
            label.set_horizontalalignment('right')

        ax[i,2].plot(dates[date_locs], np.gradient(z[date_locs]), marker='o', ls='', color=colors[2], alpha=0.3)
        smoothed = scipy.ndimage.filters.uniform_filter(np.gradient(z[date_locs]),size=5)
        ax[i,2].plot(dates[date_locs], smoothed, lw=2, color=colors[2], zorder=0)
        for label in ax[i,2].get_xticklabels():
            label.set_rotation(30)
            label.set_horizontalalignment('right')

            
        ax[i,0].grid(True, alpha=0.4)
        ax[i,1].grid(True, alpha=0.4)
        ax[i,2].grid(True, alpha=0.4)

        if i == 0:
            ax[i,0].set_title('dx/dt', color=colors[0], fontsize=16)
            ax[i,1].set_title('dy/dt', color=colors[1], fontsize=16)
            ax[i,2].set_title('dz/dt', color=colors[2], fontsize=16)
    plt.tight_layout()
    pdb.set_trace()

if __name__ == '__main__':
    target = '2MASS J23255604-0259508'
    drift_rate(target)