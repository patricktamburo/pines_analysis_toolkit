from pines_analysis_toolkit.utils import pines_dir_check, short_name_creator, pines_log_reader, get_source_names

from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from scipy.stats import pearsonr, sigmaclip
from scipy.interpolate import CubicSpline


import matplotlib.pyplot as plt 
colormap = plt.cm.viridis
import matplotlib.dates as mdates

from sklearn import linear_model
import numpy as np 
import julian 
import os
from natsort import natsorted
import pandas as pd
import itertools as it
from glob import glob
import time 
import julian 
import shutil
from copy import deepcopy
from pathlib import Path
import fileinput

def standard_x_range(times):
    """Calculates a standard x range so that light curves taken on different nights are shown on the same scale.

    :param times: Array of time stamps
    :type times: numpy array
    :return: Maximum time span of the nights in the times array in units of days
    :rtype: float
    """
    time_durations = np.zeros(len(times))
    for i in range(len(times)):
        time_durations[i] = (times[i][-1] - times[i][0])*24 + 0.25 #Gets the span of the night's observations in **HOURS**. Adds on an additional .25 hours for padding.

    return np.ceil(np.max(time_durations))/24 #Returns maximum time span of the nights in **DAYS** 

def standard_y_range(y, multiple=3.5):
    """Calculates a standard y range so that outliers don't affect the ylims of light curves.

    :param y: Array of flux measurements
    :type y: numpy array
    :param multiple: number of standard deviations to use as the y range, defaults to 3.5
    :type multiple: float, optional
    :return: A y range for the data. You should plot with, i.e., plt.ylim(1-the_return_value, 1+the_return_value)
    :rtype: float
    """
    stds = np.zeros(len(y))
    for i in range(len(y)):
        stds[i] = sigma_clipped_stats(y[i])[2]
    return multiple*np.nanmax(stds)

def raw_flux_plot(phot_path, mode='night'): 
    """Creates a plot of raw flux versus time for the target and reference stars. 
    Plot will save to same diretory as phot_path.

    :param phot_path: Path to the photometry csv file. 
    :type phot_path: pathlib.PosixPath
    :param mode: whether to run in 'night' or 'global' mode, defaults to 'night'
    :type mode: str, optional
    """

    analysis_path = phot_path.parent.parent/('analysis')
    if 'aper' in phot_path.name:
        phot_type = 'aper'
    else:
        raise RuntimeError('Have not implemented psf photometry yet!')

    if mode != 'night' and mode != 'global':
        raise ValueError("mode must be either 'night' or 'global'.")

    ap_rad = phot_path.name.split('_')[4]+'_'+phot_path.name.split('_')[1]

    phot_df = pines_log_reader(phot_path)
    times = np.array(phot_df['Time BJD TDB'])
    if mode == 'night':
        night_inds = night_splitter(times)
    #Or treat all observations as a single "night".
    elif mode == 'global':
        night_inds = [list(np.arange(len(times)))]
    num_nights = len(night_inds)

    sources = get_source_names(phot_df)
    refs = sources[1:]
    num_refs = len(refs)

    cm = plt.get_cmap('viridis_r')
    markers = ['+', 'x', '*', 'X']

    broken_times = []
    for i in range(num_nights):
        broken_times.append(np.array(times[night_inds[i]]))
    standard_time_range = standard_x_range(broken_times)

    fig, axis = plt.subplots(nrows=1, ncols=num_nights, figsize=(17,5), sharey=True)
    
    for i in range(num_nights):
        inds = night_inds[i]
        raw_targ_flux = np.array(phot_df[sources[0]+' Flux'][inds], dtype='float')
        raw_targ_err = np.array(phot_df[sources[0]+' Flux Error'][inds], dtype='float')
        raw_ref_flux = np.zeros((len(inds), num_refs))
        raw_ref_err  = np.zeros((len(inds), num_refs))
        for j in range(num_refs):
            n = refs[j]
            raw_ref_flux[:,j] = phot_df[n+' Flux'][inds]
            raw_ref_err[:,j] = phot_df[n+' Flux Error'][inds]

        #Set the axis behavior depending if there are 1 or more nights. 
        if num_nights == 1:
            ax = axis
        elif num_nights > 1:
            ax = axis[i]

        ax.set_yscale('log')

        #Plot the references. 
        for j in range(num_refs):
            color = cm(1.*j/num_refs) 
            marker = markers[j % len(markers)]
            ax.errorbar(times[inds], raw_ref_flux[:,j], raw_ref_err[:,j], marker=marker, linestyle='', label='Ref. '+str(j+1), color=color)
        
        #Set labels. 
        if i == 0:
            ax.set_ylabel('Photons', fontsize=20)
        ax.tick_params(labelsize=16)
        ax.set_xlabel('Time (BJD$_{TDB}$)', fontsize=20)

        #Plot the target. 
        ax.errorbar(times[inds], raw_targ_flux, raw_targ_err, marker='o', linestyle='-', label='Target', mfc='none', mew=2, color='tab:orange')

        if i == num_nights - 1:
            ax.legend(bbox_to_anchor=(1.03, 1.1), fontsize=10)
        
        ax.set_xlim(np.mean(times[inds])-standard_time_range/2, np.mean(times[inds])+standard_time_range/2)
        ax.grid(alpha=0.2)

    plt.suptitle(sources[0]+' Nightly Raw Flux', fontsize=20)
    plt.subplots_adjust(left=0.07, wspace=0.05, top=0.92, bottom=0.17)

    if phot_type == 'aper':
        filename = mode+'_raw_flux.png'
        if not os.path.exists(analysis_path/('aper_phot_analysis/'+ap_rad)):
            os.mkdir(analysis_path/('aper_phot_analysis/'+ap_rad))
        output_path = analysis_path/('aper_phot_analysis/'+ap_rad+'/'+filename)
        print('Saving {} raw flux plot to {}.'.format(mode, output_path))
        plt.savefig(output_path, dpi=300)
    else:
        raise RuntimeError("HAVE TO UPDATE!!!")

def alc_plot(weighted_lc_path, mode='night'):
    """Creates a standard plot of normalized target and ALC flux. 
    Plot will save to same directory as weighted_lc_path.

    :param weighted_lc_path: path to a weighted lc csv (output from weighted_lightcurve)
    :type weighted_lc_path: pathlib.PosixPath
    :param mode: whether to run in 'night' or 'global' mode, defaults to 'night'
    :type mode: str, optional
    :param force_output_path: [description], defaults to ''
    :type force_output_path: str, optional
    """
        
    weighted_lc_df = pines_log_reader(weighted_lc_path)
    sources = get_source_names(weighted_lc_df)

    times = np.array(weighted_lc_df['Time BJD TDB'])

    if mode == 'night':
        night_inds = night_splitter(times)
    #Or treat all observations as a single "night".
    elif mode == 'global':
        night_inds = [list(np.arange(len(times)))]
    num_nights = len(night_inds)

    norm_targ_flux = np.array(weighted_lc_df[sources[0]+' Normalized Flux'], dtype='float')
    norm_targ_err = np.array(weighted_lc_df[sources[0]+' Normalized Flux Error'], dtype='float')
    alc_flux = np.array(weighted_lc_df['ALC Flux'], dtype='float')
    alc_err = np.array(weighted_lc_df['ALC Flux Error'], dtype='float')
    broken_times = []
    broken_flux = []
    for i in range(num_nights):
        broken_times.append(np.array(times[night_inds[i]]))
        broken_flux.append(np.array(norm_targ_flux[night_inds[i]]))
    standard_time_range = standard_x_range(broken_times)
    standard_y = standard_y_range(broken_flux, multiple=4.0)

    fig, axis = plt.subplots(nrows=1, ncols=num_nights, figsize=(17,5), sharey=True)
    plt.subplots_adjust(left=0.07, wspace=0.05, top=0.92, bottom=0.17)

    for i in range(num_nights):
        inds = night_inds[i]
        #Set the axis behavior depending if there are 1 or more nights. 
        if num_nights == 1:
            ax = axis
        elif num_nights > 1:
            ax = axis[i]
        
        if i == 0:
            ax.set_ylabel('Normalized Flux', fontsize=20)
        ax.tick_params(labelsize=16)
        ax.set_xlabel('Time (BJD$_{TDB}$)', fontsize=20)

        if np.nanstd(norm_targ_flux[inds]) > np.nanstd(alc_flux[inds]):
            targ_zorder = 1
            alc_zorder = 2
        else:
            targ_zorder= 2
            alc_zorder= 1
        
        #Plot the target. 
        ax.errorbar(times[inds], norm_targ_flux[inds], norm_targ_err[inds], color='tab:orange', mfc='none', marker='.', ms=7, lw=1.5, linestyle='', label='Target', mew=2, zorder=targ_zorder)
        
        #Plot the ALC. 
        ax.errorbar(times[inds], alc_flux[inds], alc_err[inds], color='tab:blue', lw=2, marker='.', linestyle='', label='ALC', zorder=alc_zorder, mfc='none', mew=2, ms=5)
        
        if i == num_nights - 1:
            ax.legend(fontsize=12)

        ax.set_xlim(np.mean(times[inds])-standard_time_range/2, np.mean(times[inds])+standard_time_range/2)
        ax.set_ylim(1-standard_y, 1+standard_y)
        ax.grid(alpha=0.2)
    plt.suptitle(sources[0]+' Nightly Normalized Flux', fontsize=20)

    #Save the figure. 
    output_path = weighted_lc_path.parent/(mode+'_normalized_flux.png')
    print('Saving {} normalized flux plot to {}.'.format(mode, output_path))
    plt.savefig(output_path, dpi=300)

def corr_target_plot(weighted_lc_path, mode='night'):
    """Creates a plot of the target light curve using output from weighted_lightcurve.

    :param weighted_lc_path: path to weighted light curve csv
    :type weighted_lc_path: pathlib.PosixPath
    :param mode: whether to run in 'night' or 'global' normalization mode, defaults to 'night'
    :type mode: str, optional
    """
    weighted_lc_df = pines_log_reader(weighted_lc_path)
    sources = get_source_names(weighted_lc_df)

    times = np.array(weighted_lc_df['Time BJD TDB'])

    if mode == 'night':
        night_inds = night_splitter(times)
    #Or treat all observations as a single "night".
    elif mode == 'global':
        night_inds = [list(np.arange(len(times)))]
    num_nights = len(night_inds)

    corr_targ_flux = np.array(weighted_lc_df[sources[0]+' Corrected Flux'], dtype='float')
    broken_times = []
    broken_flux = []
    for i in range(num_nights):
        broken_times.append(np.array(times[night_inds[i]]))
        broken_flux.append(np.array(corr_targ_flux[night_inds[i]]))
    standard_time_range = standard_x_range(broken_times)
    standard_y = standard_y_range(broken_flux, multiple=3.0)
    
    #Plot the corrected target flux. 
    fig, axis = plt.subplots(nrows=1,ncols=num_nights,figsize=(17,5), sharey=True)
    #plt.subplots_adjust(left=0.07, wspace=0.05, top=0.86, bottom=0.17, right=0.96)
    plt.subplots_adjust(left=0.07, wspace=0.05, top=0.92, bottom=0.17)

    for i in range(num_nights):
        inds = night_inds[i]
        ut_date = julian.from_jd(times[inds][0])
        ut_date_str = 'UT '+ut_date.strftime('%b. %d %Y')

        #Set the axis behavior depending if there are 1 or more nights. 
        if num_nights == 1:
            ax = axis
        elif num_nights > 1:
            ax = axis[i]

        binned_times, binned_flux, binned_errs = block_binner(times[inds], corr_targ_flux[inds])        
        ax.axhline(1, color='r', lw=2, zorder=0)
        ax.grid(alpha=0.2)
        #Plot the corrected target and corrected binned target flux. 
        ax.plot(times[inds], corr_targ_flux[inds], linestyle='', marker='o', zorder=1, color='darkgrey', ms=5)
        ax.errorbar(binned_times, binned_flux, binned_errs, linestyle='', marker='o', zorder=2, color='k', capsize=2, ms=7)
        
        #Plot the 5-sigma decrement line. 
        five_sig = 5*np.nanmean(binned_errs)
        ax.axhline(1-five_sig, color='k', lw=2, linestyle='--', zorder=0)
        ax.tick_params(labelsize=16)
        
        #Set labels.
        if i == 0:
            ax.set_ylabel('Normalized Flux', fontsize=20)
        ax.set_xlabel('Time (BJD$_{TDB}$)', fontsize=20)
        
        ax.set_xlim(np.mean(times[inds])-standard_time_range/2, np.mean(times[inds])+standard_time_range/2)
        
        ax.set_ylim(1-standard_y, 1+standard_y)
        if mode == 'night':
            ax.text(np.mean(times[inds]), ax.get_ylim()[1],ut_date_str, fontsize=18,ha='center',va='top')

    plt.suptitle(sources[0], fontsize=20)

    #Save the figure. 
    output_path = weighted_lc_path.parent/(mode+'_weighted_lc.png')
    print('Saving {} normalized flux plot to {}.'.format(mode, output_path))
    plt.savefig(output_path, dpi=300)

def corr_all_sources_plot(weighted_lc_path, bin_mins=0.0, force_y_range=0.1):
    """Creates corrected light curves for all sources. Saves to 'corr_ref_plots' directory in the same directory as weighted_lc_path.

    :param weighted_lc_path: path to weighted light curve csv, output from weighted_lightcurve
    :type weighted_lc_path: pathlib.PosixPath
    :param bin_mins: number of minutes to bin data over for staring observations, defaults to 0.0
    :type bin_mins: float, optional
    :param force_y_range: force y range for all light curves for easier comparison, defaults to 0.1
    :type force_y_range: float, optional
    """
    print('Generating corrected flux plots for all sources...\n')

    data = pines_log_reader(weighted_lc_path)
    output_path = weighted_lc_path.parent/('corr_ref_plots/')
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    sources = get_source_names(data)
    ref_names = sources[1:]
    num_refs = len(ref_names)
    
    times = np.array(data['Time BJD TDB'])
    night_inds = night_splitter(times)
    num_nights = len(night_inds)
    
    cmap = plt.get_cmap('viridis')
    for i in range(num_refs + 1):
        fig, ax = plt.subplots(nrows=1, ncols=num_nights, figsize=(17,5), sharey=True)
        plt.subplots_adjust(left=0.07, wspace=0.05, top=0.92, bottom=0.17)

        if i == 0:
            color = cmap(0)
            flux = np.array(data[sources[0]+' Corrected Flux'], dtype='float64')
            flux_err = np.array(data[sources[0]+' Corrected Flux Error'], dtype='float64')
            title = sources[0]
            output_name = sources[0]+'_corrected_flux.png'

        else:
            color = cmap(95)
            ref_name = ref_names[i-1]
            flux = np.array(data[ref_name+' Corrected Flux'], dtype='float64')
            flux_err = np.array(data[ref_name+' Corrected Flux Error'], dtype='float64')
            num = ref_names[i-1].split(' ')[1].zfill(2)
            output_name = 'reference_'+num+'_corrected_flux.png'

        for j in range(num_nights):
            if i != 0: 
                weight = np.array(data[ref_name+' ALC Weight'])[night_inds[j]][0]
                title = ref_name.replace('erence','.')+', weight = {:1.3f}'.format(weight)

            if j == 0:
                if num_nights == 1:
                    ax.set_ylabel('Normalized Flux', fontsize=20)
                else:
                    ax[j].set_ylabel('Normalized Flux', fontsize=20)
            
            inds = night_inds[j]

            binned_time, binned_flux, binned_err = block_binner(times[inds], flux[inds], bin_mins=bin_mins)          

            if num_nights == 1:
                ax.plot(times[inds], flux[inds], color=color, linestyle='', marker='.', alpha=0.25)
                ax.errorbar(binned_time, binned_flux, binned_err, color=color, linestyle='', marker='o', ms=10, mfc='none', mew=2)
                ax.set_xlabel('Time (BJD$_{TDB}$)', fontsize=20)
                ax.tick_params(labelsize=16)
                ax.axhline(1, color='k', alpha=0.7, lw=1, zorder=0)
                ax.grid(alpha=0.2)
                ax.set_title(title, fontsize=16, color=color)
                ax.set_ylim(1-force_y_range,1+force_y_range)
            else:
                ax[j].plot(times[inds], flux[inds], color=color, linestyle='', marker='.', alpha=0.25)
                ax[j].errorbar(binned_time, binned_flux, binned_err, color=color, linestyle='', marker='o', ms=10, mfc='none', mew=2)
                ax[j].set_xlabel('Time (BJD$_{TDB}$)', fontsize=20)
                ax[j].tick_params(labelsize=16)
                ax[j].axhline(1, color='k', alpha=0.7, lw=1, zorder=0)
                ax[j].grid(alpha=0.2)
                ax[j].set_title(title, fontsize=16, color=color)
                ax[j].set_ylim(1-force_y_range,1+force_y_range)

        plt.savefig(output_path/output_name, dpi=300)
        plt.close()
    
def block_binner(raw_times, raw_flux, time_threshold=0.15, bin_mins=0.0):
    """Bins PINES data over blocks.

    :param raw_times: array of times
    :type raw_times: numpy array
    :param raw_flux: array of fluxes
    :type raw_flux: numpy array
    :param time_threshold: the gap between observations, in hours, above which sets of observations will be considered different blocks, defaults to 0.15
    :type time_threshold: float, optional
    :param bin_mins: minutes over which to bin data for staring observations, defaults to 0.0
    :type bin_mins: float, optional
    :return: arrays of binned times, flux, and errors
    :rtype: numpy arrays
    """
    block_inds = block_splitter(raw_times, time_threshold=time_threshold, bin_mins=bin_mins)
    n_blocks = len(block_inds)
    bin_times      = np.zeros(n_blocks)
    bin_flux       = np.zeros(n_blocks)
    bin_flux_err   = np.zeros(n_blocks)
    for i in range(n_blocks):
        inds = block_inds[i]
        bin_times[i]    = np.nanmean(raw_times[inds])
        bin_flux[i]     = np.nanmean(raw_flux[inds])
        bin_flux_err[i] = np.nanstd(raw_flux[inds])/np.sqrt(len(inds))
    return bin_times, bin_flux, bin_flux_err

def block_splitter(times, time_threshold=0.15, bin_mins=0.0):
    """Finds individual blocks of data given a single night of exposture times. 

    :param times: array of times
    :type times: numpy array
    :param time_threshold: the gap between observations, in hours, above which sets of observations will be considered different blocks, defaults to 0.15
    :type time_threshold: float, optional
    :param bin_mins: can alternatively choose a time duration over which to bin data, which is needed for staring data, defaults to 0.0. 
    :type bin_mins: float, optional
    :return: list of length n_blocks, each entry containing the indices of data for that block
    :rtype: list
    """
    times = np.array(times) 

    #Staring observations. TODO: This will not work if there is a mix of staring/hopping observations on a single night!
    if bin_mins != 0.0:
        time_bin = bin_mins #Minutes over which to bin. 
        block_boundaries = np.where(np.gradient(((times - times[0]) * (24*60)) % time_bin) < 0.2)[0]
    else:
        block_boundaries = np.where(np.gradient(times) > time_threshold/24)[0]

    num_blocks = int(1 + len(block_boundaries) / 2)
    block_inds = [[] for x in range(num_blocks)]
    for j in range(num_blocks):
        if j == 0:
            block_inds[j].extend(np.arange(0,block_boundaries[0]+1))
        elif j > 0 and j < num_blocks - 1:
            block_inds[j].extend(np.arange(block_boundaries[2*j-1],block_boundaries[2*j]+1))
        else:
            block_inds[j].extend(np.arange(block_boundaries[2*j-1],len(times)))
    
    return block_inds

def plot_style():
    """Sets some font sizes for standard plots. 

    :return: plot style parameters
    :rtype: floats
    """

    title_size = 20 
    axis_title_size = 20
    axis_ticks_font_size = 16
    legend_font_size = 12
    return (title_size, axis_title_size, axis_ticks_font_size, legend_font_size)

def seeing_plot(short_name, bin_mins=0.0, force_output_path=''):
    """Creates a plot of seeing versus time.

    :param short_name: short name of the target
    :type short_name: str
    :param bin_mins: number of minutes to bin over data for staring observations, defaults to 0.0
    :type bin_mins: float, optional
    :param force_output_path: user-chosen path if you do not want to use the default ~/Documents/PINES_analysis_toolkit/ directory for analysis, defaults to ''
    :type force_output_path: str, optional
    """
    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()

    #Get plot style parameters. 
    title_size, axis_title_size, axis_ticks_font_size, legend_font_size = plot_style()
    
    centroided_sources_path = pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv')
    if not os.path.exists(centroided_sources_path):
        print('ERROR: no centroided sources csv, not doing seeing plot.')
    else:
        centroided_sources = pines_log_reader(centroided_sources_path)

    #Get list of souce names in the centroid output.
    source_names = get_source_names(centroided_sources)
    centroided_sources.columns = centroided_sources.keys().str.strip()

    #Read in times and seeing values. 
    times_full = np.array(centroided_sources['Time BJD TDB'], dtype='float')
    seeing = np.array(centroided_sources['Seeing'], dtype='float')

    #Split up times by nights.
    night_inds = night_splitter(times_full)
    num_nights = len(night_inds)
    times_nights = [times_full[night_inds[i]] for i in range(num_nights)]
    standard_x = standard_x_range(times_nights)

    fig, ax = plt.subplots(nrows=1, ncols=num_nights, figsize=(17,5), sharey=True)
    for i in range(num_nights):
        if i == 0:
            if num_nights == 1:
                ax.set_ylabel('Seeing (")', fontsize=axis_title_size)
            else:
                ax[i].set_ylabel('Seeing (")', fontsize=axis_title_size)

        inds = night_inds[i]
        

        #bin
        block_inds = block_splitter(times_nights[i], bin_mins=bin_mins)
        block_x = np.zeros(len(block_inds))
        block_y = np.zeros(len(block_inds))
        block_y_err = np.zeros(len(block_inds))
        for j in range(len(block_inds)):
            block_x[j] = np.nanmean(times_nights[i][block_inds[j]])
            block_y[j] = np.nanmean(seeing[inds][block_inds[j]])
            block_y_err[j] = np.nanstd(seeing[inds][block_inds[j]]) / np.sqrt(len(seeing[inds][block_inds[j]]))
       
        
        #Interpolate each night's seeing.
        fit_times = np.linspace(block_x[0], block_x[-1], 1000)
        interp = CubicSpline(block_x, block_y)
        interp_fit = interp(fit_times)
        
        if num_nights == 1:
            ax.plot(times_nights[i], seeing[inds], marker='.', linestyle='', alpha=0.3, label='Raw seeing')
            ax.tick_params(labelsize=axis_ticks_font_size)
            ax.set_xlabel('Time (BJD$_{TDB}$)', fontsize=axis_title_size)
            ax.grid(alpha=0.2)
            ax.set_xlim(np.mean(times_nights[i])-standard_x/2, np.mean(times_nights[i])+standard_x/2)
            ax.errorbar(block_x, block_y, block_y_err, marker='o', linestyle='', color='tab:blue', ms=8, mfc='none', mew=2, label='Bin seeing')
            ax.plot(fit_times, interp_fit, color='r', lw=2, zorder=0, alpha=0.7, label='CS Interp.')
        else:
            ax[i].plot(times_nights[i], seeing[inds], marker='.', linestyle='', alpha=0.3, label='Raw seeing')
            ax[i].tick_params(labelsize=axis_ticks_font_size)
            ax[i].set_xlabel('Time (BJD$_{TDB}$)', fontsize=axis_title_size)
            ax[i].grid(alpha=0.2)
            ax[i].set_xlim(np.mean(times_nights[i])-standard_x/2, np.mean(times_nights[i])+standard_x/2)
            ax[i].errorbar(block_x, block_y, block_y_err, marker='o', linestyle='', color='tab:blue', ms=8, mfc='none', mew=2, label='Bin seeing')
            ax[i].plot(fit_times, interp_fit, color='r', lw=2, zorder=0, alpha=0.7, label='CS Interp.')
    if num_nights == 1:
        ax.legend(bbox_to_anchor=(1.01, 0.5), fontsize=legend_font_size)
    else:
        ax[i].legend(bbox_to_anchor=(1.01, 0.5), fontsize=legend_font_size)
    plt.suptitle(short_name+' Seeing Measurements', fontsize=title_size)
    plt.subplots_adjust(left=0.07, wspace=0.05, top=0.92, bottom=0.17)

    output_filename = pines_path/('Objects/'+short_name+'/analysis/diagnostic_plots/'+short_name+'_seeing.png')
    print('Saving seeing plot to {}.'.format(output_filename))
    plt.savefig(output_filename, dpi=300)
    return

def relative_cutout_position_plot(short_name, force_output_path=''):
    """Creates a plot of relative cutout positions versus time for all sources.

    :param short_name: short name for the target
    :type short_name: str
    :param force_output_path: user-chosen path if you do not want to use the default ~/Documents/PINES_analysis_toolkit/ directory for analysis, defaults to ''
    :type force_output_path: str, optional
    """
    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()
    
    centroided_sources_path = pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv')
    if not os.path.exists(centroided_sources_path):
        print('ERROR: no centroided sources csv, not doing relative cutout position plot.')
    else:
        centroided_sources = pines_log_reader(centroided_sources_path)

    #Get plot style parameters. 
    title_size, axis_title_size, axis_ticks_font_size, legend_font_size = plot_style()
    
    #Get list of souce names in the centroid output.
    source_names = get_source_names(centroided_sources)
    centroided_sources.columns = centroided_sources.keys().str.strip()

    #Get times from the centroid output and split them by night. 
    times_full = np.array(centroided_sources['Time BJD TDB'], dtype='float')
    night_inds = night_splitter(times_full)
    num_nights = len(night_inds)
    times_nights = [times_full[night_inds[i]] for i in range(num_nights)]
    standard_x = standard_x_range(times_nights)

    #Get the box size (I don't like that this is being determined by using the mean of the data...output it from centroider?)    
    box_w = int(np.round(2*np.nanmean(np.array(centroided_sources['Reference 1 Cutout X'], dtype='float')),0))
        
    fig, ax = plt.subplots(nrows=2, ncols=num_nights, figsize=(17,9), sharey=True)
    plt.subplots_adjust(left=0.07, hspace=0.05, wspace=0.05, top=0.92, bottom=0.17)
    markers = ['+', 'x', '*', 'X']
    for j in range(num_nights):
        inds = night_inds[j]
        if j == 0: 
            if num_nights == 1:
                ax[0].set_ylabel('Cutout X Position', fontsize=axis_title_size)
                ax[1].set_ylabel('Cutout Y Position', fontsize=axis_title_size)
                #ax[0].set_ylim(0.25,0.75)
            else:
                ax[0,j].set_ylabel('Cutout X Position', fontsize=axis_title_size)
                ax[1,j].set_ylabel('Cutout Y Position', fontsize=axis_title_size)
                #ax[0,j].set_ylim(0.25,0.75)
        for i in range(len(source_names)):
            if i == 1:
                continue
            cutout_x = np.array(centroided_sources[source_names[i]+' Cutout X'][inds], dtype='float')
            cutout_y = np.array(centroided_sources[source_names[i]+' Cutout Y'][inds], dtype='float')
            if i == 0:
                marker = 'o'
                label = 'Target'
            else:
                marker = markers[(i-1) % len(markers)]
                label = 'Ref. '+str(i)

            if num_nights == 1:
                ax[0].plot(times_nights[j], cutout_x, marker=marker, label=label, linestyle='')
                ax[1].plot(times_nights[j], cutout_y, marker=marker, linestyle='')
            else:
                ax[0,j].plot(times_nights[j], cutout_x, marker=marker, label=label, linestyle='')
                ax[1,j].plot(times_nights[j], cutout_y, marker=marker, linestyle='')

        if num_nights == 1:
            ax[0].tick_params(labelsize=axis_ticks_font_size)
            ax[0].set_xticklabels([])
            ax[0].axhline(box_w/2, zorder=0, color='r', label='Center pix.', lw=2)
            ax[0].set_xlim(np.mean(times_nights[j])-standard_x/2, np.mean(times_nights[j])+standard_x/2)
            ax[0].grid(alpha=0.2)
            ax[1].tick_params(labelsize=axis_ticks_font_size)
            ax[1].axhline(box_w/2, zorder=0, color='r', label='Center pix.', lw=2)
            ax[1].set_xlim(np.mean(times_nights[j])-standard_x/2, np.mean(times_nights[j])+standard_x/2)
            ax[1].set_xlabel('Time (BJD$_{TDB}$)', fontsize=axis_title_size)
            ax[1].grid(alpha=0.2)
        else:
            ax[0,j].tick_params(labelsize=axis_ticks_font_size)
            ax[0,j].set_xticklabels([])
            ax[0,j].axhline(box_w/2, zorder=0, color='r', label='Center pix.', lw=2)
            ax[0,j].set_xlim(np.mean(times_nights[j])-standard_x/2, np.mean(times_nights[j])+standard_x/2)
            ax[0,j].grid(alpha=0.2)
            ax[1,j].tick_params(labelsize=axis_ticks_font_size)
            ax[1,j].axhline(box_w/2, zorder=0, color='r', label='Center pix.', lw=2)
            ax[1,j].set_xlim(np.mean(times_nights[j])-standard_x/2, np.mean(times_nights[j])+standard_x/2)
            ax[1,j].set_xlabel('Time (BJD$_{TDB}$)', fontsize=axis_title_size)
            ax[1,j].grid(alpha=0.2)

        if j == num_nights - 1:
            if num_nights == 1:
                ax[0].legend(bbox_to_anchor=(1.01, 1.0), fontsize=legend_font_size)
            else:
                ax[0,j].legend(bbox_to_anchor=(1.01, 1.0), fontsize=legend_font_size)

    plt.suptitle(short_name+' Cutout Centroid Positions', fontsize=title_size)

    output_filename = pines_path/('Objects/'+short_name+'/analysis/diagnostic_plots/'+short_name+'_cutout_positions.png')
    print('Saving relative cutout position plot to {}.'.format(output_filename))
    plt.savefig(output_filename, dpi=300)

    return

def absolute_image_position_plot(short_name, bin_mins=0.0, force_output_path=''):
    """Creates a plot of target x/y pixel positions versus time. 

    :param short_name: short name of the target
    :type short_name: str
    :param bin_mins:  number of minutes to bin over data for staring observations, defaults to 0.0
    :type bin_mins: float, optional
    :param force_output_path: user-chosen path if you do not want to use the default ~/Documents/PINES_analysis_toolkit/ directory for analysis, defaults to ''
    :type force_output_path: str, optional
    """
    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()

    centroided_sources_path = pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv')
    if not os.path.exists(centroided_sources_path):
        print('ERROR: no centroided sources csv, not doing absolute position plot.')
    else:
        centroided_sources = pines_log_reader(centroided_sources_path)

    #Get plot style parameters. 
    title_size, axis_title_size, axis_ticks_font_size, legend_font_size = plot_style()
    
    #Get list of souce names in the centroid output.
    source_names = get_source_names(centroided_sources)
    centroided_sources.columns = centroided_sources.keys().str.strip()

    #Get times from the centroid output and split them by night. 
    times_full = np.array(centroided_sources['Time BJD TDB'])
    night_inds = night_splitter(times_full)
    num_nights = len(night_inds)
    times_nights = [times_full[night_inds[i]] for i in range(num_nights)]
    standard_x = standard_x_range(times_nights)
    
    targ_ind = np.where(['2MASS' in i for i in source_names])[0][0]
    targ_name = source_names[targ_ind]
    fig, ax = plt.subplots(nrows=2, ncols=num_nights, figsize=(17,9), sharex='col', sharey='row')
    plt.subplots_adjust(left=0.07, hspace=0.05, wspace=0.05, top=0.92, bottom=0.17)
    for j in range(num_nights):
        if j == 0:
            if num_nights == 1:
                ax[0].set_ylabel('Image X', fontsize=axis_title_size)
                ax[1].set_ylabel('Image Y', fontsize=axis_title_size)
            else:
                ax[0,j].set_ylabel('Image X', fontsize=axis_title_size)
                ax[1,j].set_ylabel('Image Y', fontsize=axis_title_size)

        inds = night_inds[j]
        times = times_nights[j]
        absolute_x = np.array(centroided_sources[targ_name+' Image X'][inds], dtype='float')
        absolute_y = np.array(centroided_sources[targ_name+' Image Y'][inds], dtype='float')
        if num_nights == 1:
            ax[0].plot(times, absolute_x, marker='.', linestyle='', alpha=0.3, color='tab:blue', label='Raw x')
            ax[1].plot(times, absolute_y, marker='.', linestyle='', alpha=0.3, color='tab:orange', label='Raw y')
        else:
            ax[0,j].plot(times, absolute_x, marker='.', linestyle='', alpha=0.3, color='tab:blue', label='Raw x')
            ax[1,j].plot(times, absolute_y, marker='.', linestyle='', alpha=0.3, color='tab:orange', label='Raw y')

        #bin
        block_inds = block_splitter(times, bin_mins=bin_mins)
        block_times= np.zeros(len(block_inds))
        block_x = np.zeros(len(block_inds))
        block_x_err = np.zeros(len(block_inds))
        block_y = np.zeros(len(block_inds))
        block_y_err = np.zeros(len(block_inds))
        for k in range(len(block_inds)):
            try:
                block_times[k] = np.nanmean(times[block_inds[k]])
                block_x[k] = np.nanmean(absolute_x[block_inds[k]])
                block_x_err[k] = np.nanstd(absolute_x[block_inds[k]]) / np.sqrt(len(absolute_x[block_inds[k]]))
                block_y[k] = np.nanmean(absolute_y[block_inds[k]])
                block_y_err[k] = np.nanstd(absolute_y[block_inds[k]]) / np.sqrt(len(absolute_y[block_inds[k]]))
            except:
                pdb.set_trace()


        if num_nights == 1:
            ax[0].errorbar(block_times, block_x, block_x_err, marker='o', linestyle='', color='tab:blue', ms=8, mfc='none', mew=2, label='Bin x')
            ax[1].errorbar(block_times, block_y, block_y_err, marker='o', linestyle='', color='tab:orange', ms=8, mfc='none', mew=2, label='Bin y')
            ax[0].tick_params(labelsize=axis_ticks_font_size)
            ax[1].tick_params(labelsize=axis_ticks_font_size)
            ax[0].grid(alpha=0.2)
            ax[1].grid(alpha=0.2)
            ax[1].set_xlabel('Time (BJD$_{TDB}$)', fontsize=axis_title_size)
        else:
            ax[0,j].errorbar(block_times, block_x, block_x_err, marker='o', linestyle='', color='tab:blue', ms=8, mfc='none', mew=2, label='Bin x')
            ax[1,j].errorbar(block_times, block_y, block_y_err, marker='o', linestyle='', color='tab:orange', ms=8, mfc='none', mew=2, label='Bin y')
            ax[0,j].tick_params(labelsize=axis_ticks_font_size)
            ax[1,j].tick_params(labelsize=axis_ticks_font_size)
            ax[0,j].grid(alpha=0.2)
            ax[1,j].grid(alpha=0.2)
            ax[1,j].set_xlabel('Time (BJD$_{TDB}$)', fontsize=axis_title_size)

        if j == num_nights - 1:
            if num_nights == 1:
                ax[0].legend(bbox_to_anchor=(1.29, 1), fontsize=legend_font_size)
                ax[1].legend(bbox_to_anchor=(1.29, 1), fontsize=legend_font_size)
            else:
                ax[0,j].legend(bbox_to_anchor=(1.29, 1), fontsize=legend_font_size)
                ax[1,j].legend(bbox_to_anchor=(1.29, 1), fontsize=legend_font_size)

    plt.suptitle(targ_name+' Image Centroid Positions', fontsize=title_size)
    #plt.subplots_adjust(left=0.07, hspace=0.05, wspace=0.05, top=0.92, bottom=0.08, right=0.85)

    
   #ax.legend(bbox_to_anchor=(1.01, 1), fontsize=14)

    if num_nights == 1:
        ax[0].set_xlim(np.mean(times)-standard_x/2, np.mean(times)+standard_x/2)
        ax[1].set_xlim(np.mean(times)-standard_x/2, np.mean(times)+standard_x/2)
    else:
        ax[0,j].set_xlim(np.mean(times)-standard_x/2, np.mean(times)+standard_x/2)
        ax[1,j].set_xlim(np.mean(times)-standard_x/2, np.mean(times)+standard_x/2)

    output_filename = pines_path/('Objects/'+short_name+'/analysis/diagnostic_plots/'+targ_name+'_image_positions.png')
    print('Saving absolute image position plot to {}.'.format(output_filename))
    plt.savefig(output_filename, dpi=300)

    return 

def background_plot(short_name, gain=8.21, bin_mins=0.0, force_output_path=''):
    """Creates a plot of measured target background versus time.

    :param short_name: short name of the target
    :type short_name: str
    :param gain: gain of the detector in e-/ADU, defaults to 8.21
    :type gain: float, optional
    :param bin_mins: number of minutes to bin over data for staring observations, defaults to ''
    :type bin_mins: str, optional
    :param force_output_path: user-chosen path if you do not want to use the default ~/Documents/PINES_analysis_toolkit/ directory for analysis, defaults to ''
    :type force_output_path: str, optional
    """
    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()


    #Get plot style parameters. 
    title_size, axis_title_size, axis_ticks_font_size, legend_font_size = plot_style()

    analysis_path = pines_path/('Objects/'+short_name+'/analysis')
    phot_path = pines_path/('Objects/'+short_name+'/aper_phot')
    phot_files = np.array(natsorted([x for x in phot_path.glob('*.csv')]))

    if os.path.exists(analysis_path/('optimal_aperture.txt')):
        with open(analysis_path/('optimal_aperture.txt'), 'r') as f:
            best_ap = f.readlines()[0].split(':  ')[1].split('_')[0]
        ap_list = np.array([str(i).split('/')[-1].split('_')[4] for i in phot_files])
        best_ap_ind = np.where(ap_list == best_ap)[0][0]
    else:
        print('No optimal_aperture.txt file for {}.\nUsing first photometry file in {}.'.format(target, phot_path))
        best_ap_ind = 0

    phot_file = phot_files[best_ap_ind]
    phot_df = pines_log_reader(phot_file)

    backgrounds = np.array(phot_df[short_name+' Background'], dtype='float')/gain
    times_full = np.array(phot_df['Time BJD TDB'], dtype='float')
    night_inds = night_splitter(times_full)
    num_nights = len(night_inds)
    times_nights = [times_full[night_inds[i]] for i in range(num_nights)]
    standard_x = standard_x_range(times_nights)

    fig, ax = plt.subplots(nrows=1, ncols=num_nights, figsize=(17,5), sharey=True)
    for i in range(num_nights):
        if i == 0:
            if num_nights == 1:
                ax.set_ylabel('Background (ADU)', fontsize=axis_title_size)
            else:
                ax[i].set_ylabel('Background (ADU)', fontsize=axis_title_size)

        inds = night_inds[i]
        #bin
        block_inds = block_splitter(times_full[inds], bin_mins=bin_mins)
        block_x = np.zeros(len(block_inds))
        block_y = np.zeros(len(block_inds))
        block_y_err = np.zeros(len(block_inds))
        for j in range(len(block_inds)):
            block_x[j] = np.nanmean(times_full[inds][block_inds[j]])
            block_y[j] = np.nanmean(backgrounds[inds][block_inds[j]])
            block_y_err[j] = np.nanstd(backgrounds[inds][block_inds[j]]) / np.sqrt(len(backgrounds[inds][block_inds[j]]))

        block_x = block_x[~np.isnan(block_y)]
        block_y_err = block_y_err[~np.isnan(block_y)]
        block_y = block_y[~np.isnan(block_y)]
        
        #Interpolate each night's seeing.
        fit_times = np.linspace(block_x[0], block_x[-1], 1000)
        try:
            interp = CubicSpline(block_x, block_y)
        except:
            pdb.set_trace()
        interp_fit = interp(fit_times)

        if num_nights == 1:
            ax.plot(times_full[inds], backgrounds[inds], marker='.', linestyle='', color='tab:orange', alpha=0.3, label='Raw bkg.')
            ax.tick_params(labelsize=axis_ticks_font_size)
            ax.set_xlabel('Time (BJD$_{TDB}$)', fontsize=axis_title_size)
            ax.grid(alpha=0.2)
            ax.set_xlim(np.mean(times_full[inds])-standard_x/2, np.mean(times_full[inds]+standard_x/2))
            ax.errorbar(block_x, block_y, block_y_err, marker='o', linestyle='', color='tab:orange', ms=8, mfc='none', mew=2, label='Bin bkg.')
            ax.legend(bbox_to_anchor=(1.01, 0.5), fontsize=legend_font_size)
            ax.plot(fit_times, interp_fit, color='b', lw=2, zorder=0, alpha=0.7, label='CS Interp.')
        else:
            ax[i].plot(times_full[inds], backgrounds[inds], marker='.', linestyle='', color='tab:orange', alpha=0.3, label='Raw bkg.')
            ax[i].tick_params(labelsize=axis_ticks_font_size)
            ax[i].set_xlabel('Time (BJD$_{TDB}$)', fontsize=axis_title_size)
            ax[i].grid(alpha=0.2)
            ax[i].set_xlim(np.mean(times_full[inds])-standard_x/2, np.mean(times_full[inds]+standard_x/2))
            ax[i].errorbar(block_x, block_y, block_y_err, marker='o', linestyle='', color='tab:orange', ms=8, mfc='none', mew=2, label='Bin bkg.')
            if i == num_nights - 1:
                ax[i].legend(bbox_to_anchor=(1.01, 0.5), fontsize=legend_font_size)
            ax[i].plot(fit_times, interp_fit, color='b', lw=2, zorder=0, alpha=0.7, label='CS Interp.')

    plt.suptitle(short_name+' Background Measurements', fontsize=title_size)
    plt.subplots_adjust(left=0.07, wspace=0.05, top=0.92, bottom=0.17)
    
    output_filename = pines_path/('Objects/'+short_name+'/analysis/diagnostic_plots/'+short_name+'_backgrounds.png')
    print('Saving background plot to {}.'.format(output_filename))
    plt.savefig(output_filename, dpi=300)
    return

def airmass_plot(short_name, bin_mins=0.0, force_output_path=''):
    """Creates a plot of airmass versus time.

    :param short_name: short name of the target
    :type short_name: str
    :param bin_mins: number of minutes to bin over data for staring observations, defaults to 0.0
    :type bin_mins: float, optional
    :param force_output_path: user-chosen path if you do not want to use the default ~/Documents/PINES_analysis_toolkit/ directory for analysis, defaults to ''
    :type force_output_path: str, optional
    """
    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()
    centroided_sources_path = pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv')
    if not os.path.exists(centroided_sources_path):
        print('ERROR: no centroided sources csv, not doing airmass plot.')
    else:        
        centroided_sources = pines_log_reader(centroided_sources_path)

    #Get plot style parameters. 
    title_size, axis_title_size, axis_ticks_font_size, legend_font_size = plot_style()

    #Get list of souce names in the centroid output.
    centroided_sources.columns = centroided_sources.keys().str.strip()
    airmasses = np.array(centroided_sources['Airmass'])
    times_full = np.array(centroided_sources['Time BJD TDB'])
    night_inds = night_splitter(times_full)
    num_nights = len(night_inds)
    times_nights = [times_full[night_inds[i]] for i in range(num_nights)]
    standard_x = standard_x_range(times_nights)

    fig, ax = plt.subplots(nrows=1, ncols=num_nights, figsize=(17,5), sharey=True)
    for i in range(num_nights):
        if i == 0:
            if num_nights == 1:
                ax.set_ylabel('Airmass', fontsize=axis_title_size)
            else:
                ax[i].set_ylabel('Airmass', fontsize=axis_title_size)

        inds = night_inds[i]

        #bin
        block_inds = block_splitter(times_full[inds], bin_mins=bin_mins)
        block_x = np.zeros(len(block_inds))
        block_y = np.zeros(len(block_inds))
        block_y_err = np.zeros(len(block_inds))
        for j in range(len(block_inds)):
            block_x[j] = np.mean(times_full[inds][block_inds[j]])
            block_y[j] = np.mean(airmasses[inds][block_inds[j]])
            block_y_err[j] = np.std(airmasses[inds][block_inds[j]]) / np.sqrt(len(airmasses[inds][block_inds[j]]))

        
        
        #Interpolate each night's seeing.
        fit_times = np.linspace(block_x[0], block_x[-1], 1000)
        interp = CubicSpline(block_x, block_y)
        interp_fit = interp(fit_times)

        if num_nights == 1:
            ax.plot(times_full[inds], airmasses[inds], marker='.', linestyle='', color='m', alpha=0.3, label='Raw airmass')
            ax.tick_params(labelsize=axis_ticks_font_size)
            ax.set_xlabel('Time (BJD$_{TDB}$)', fontsize=axis_title_size)
            ax.grid(alpha=0.2)
            ax.set_xlim(np.mean(times_full[inds])-standard_x/2, np.mean(times_full[inds]+standard_x/2))
            ax.errorbar(block_x, block_y, block_y_err, marker='o', linestyle='', color='m', ms=8, mfc='none', mew=2, label='Bin airmass')
            ax.plot(fit_times, interp_fit, color='c', lw=2, zorder=0, alpha=0.7, label='CS Interp.')
        else:
            ax[i].plot(times_full[inds], airmasses[inds], marker='.', linestyle='', color='m', alpha=0.3, label='Raw airmass')
            ax[i].tick_params(labelsize=axis_ticks_font_size)
            ax[i].set_xlabel('Time (BJD$_{TDB}$)', fontsize=axis_title_size)
            ax[i].grid(alpha=0.2)
            ax[i].set_xlim(np.mean(times_full[inds])-standard_x/2, np.mean(times_full[inds]+standard_x/2))
            ax[i].errorbar(block_x, block_y, block_y_err, marker='o', linestyle='', color='m', ms=8, mfc='none', mew=2, label='Bin airmass')
            ax[i].plot(fit_times, interp_fit, color='c', lw=2, zorder=0, alpha=0.7, label='CS Interp.')
    
    if num_nights == 1:
        ax.legend(bbox_to_anchor=(1.005, 0.5), fontsize=legend_font_size)
    else:
        ax[i].legend(bbox_to_anchor=(1.005, 0.5), fontsize=legend_font_size)
    plt.suptitle(short_name+' Airmass Measurements', fontsize=title_size)
    plt.subplots_adjust(left=0.07, wspace=0.05, top=0.92, bottom=0.17)

    output_filename = pines_path/('Objects/'+short_name+'/analysis/diagnostic_plots/'+short_name+'_airmasses.png')
    print('Saving airmass plot to {}.'.format(output_filename))
    plt.savefig(output_filename, dpi=300)
    return

def night_splitter(times):
    """Finds different nights of data given a 1D array of exposure times (in days), and returns their indices. 

    :param times: array of times (in units of days)
    :type times: numpy array
    :return: list of length n_nights, each entry containing the indices of data from a particular night of observations
    :rtype: list
    """
    night_boundaries = np.where(np.gradient(times) > 7/24)[0]
    num_nights = int(1 + len(night_boundaries) / 2)
    night_inds = [[] for x in range(num_nights)]
    if num_nights == 1:
        night_inds[0].extend(np.arange(0, len(times)))
    else:
        for j in range(num_nights):
            if j == 0:
                night_inds[j].extend(np.arange(0,night_boundaries[0]+1))
            elif j > 0 and j < num_nights - 1:
                night_inds[j].extend(np.arange(night_boundaries[2*j-1],night_boundaries[2*j]+1))
            else:
                night_inds[j].extend(np.arange(night_boundaries[2*j-1],len(times)))
    return night_inds

def readnoise_calculator():
    """Reads in two bias files, subtracts them, finds the standard deviation of the difference image, then measures read noise. 
    Following http://spiff.rit.edu/classes/phys445/lectures/readout/readout.html.
    """
    G = 8.21 #e- / DN
    bias_numbers = np.arange(23,43,1) #biases from this run were images 23 - 42 on 2019/11/14. 
    unique_combos = pd.Series(list(it.combinations(np.unique(bias_numbers),2))) #Get all unique combinations of those images. 
    read_noises = np.zeros(len(unique_combos))
    for i in range(len(unique_combos)):
        bias1_path = 'P:\\Data\\201911\\20191114\\20191114.0'+str(unique_combos[i][0])+'.fits'
        bias2_path = 'P:\\Data\\201911\\20191114\\20191114.0'+str(unique_combos[i][1])+'.fits'
        bias1 = fits.open(bias1_path)[0].data
        bias2 = fits.open(bias2_path)[0].data
        diff_image = bias1-bias2 #Make the difference image. 
        read_noises[i] = np.std(diff_image)/np.sqrt(2) * G #Calculate RN, converted to e-. 
        print(read_noises[i])

def regression(flux, regressors, corr_significance=1e-2, verbose=False):
    """Does a regression of target flux against regressors, if the correlation significance is less than the significance threshold. 

    :param flux: array of flux values
    :type flux: numpy array
    :param regressors: dictionary containing labeled regressors, each a numpy array of the same length as the flux values
    :type regressors: dict
    :param corr_significance: the p-value that a correlation between the flux and an individual regressor must have to be included in the regression, defaults to 1e-2
    :type corr_significance: float, optional
    :param verbose: whether or not to print information about the regressors used, defaults to False
    :type verbose: bool, optional
    :return: regressed flux
    :rtype: numpy array
    """

    keys = np.array(list(regressors.keys()))
    good_locs = np.where(~np.isnan(flux))[0] #Only perform the regression on non-NaN values. 

    sigs = []
    if verbose:
        print('-------------------------')
        print('{:<11s} | {:>11s}'.format('Regressor', 'Signficance'))
        print('-------------------------')
    for i in range(len(regressors)):
        corr, sig = pearsonr(flux[good_locs], regressors[keys[i]][good_locs])
        sigs.append(sig)
        if verbose:
            print('{:<11s} | {:>.2e}'.format(keys[i], sig))

    use_inds = np.where(np.array(sigs) <= corr_significance)
    use_keys = keys[use_inds]

    if verbose:
        print('Using the following regressors: ')
        for i in range(len(use_keys)):
            print('{:<11s}'.format(use_keys[i]))
        print('')

    #Now set up the linear regression.
    regr = linear_model.LinearRegression()
    
    #Set up the regression dict using only the regressors with correlation significances less than corr_significance
    regress_dict = {}
    for i in range(len(use_keys)):
        regress_dict[use_keys[i]] = regressors[use_keys[i]]
    
    #Finally, add target flux
    regress_dict['flux'] = flux

    #Get list of keys
    keylist = list()
    for i in regress_dict.keys():
        keylist.append(i)

    #Create data frame of regressors.
    df = pd.DataFrame(regress_dict,columns=keylist)
    x = df[keylist[0:len(keylist)-1]]
    y = df['flux']

    if np.shape(x)[1] >0:
        regr.fit(x[~np.isnan(y)],y[~np.isnan(y)])

        #Now, define the model.
        linear_regression_model = regr.intercept_

        for i in range(len(use_keys)):
            linear_regression_model += regr.coef_[i]*regress_dict[use_keys[i]]

        #Calculate the chi^2 of the fit. 
        chi_2 = np.nansum((flux - linear_regression_model)**2/(linear_regression_model))

        #Divide out the fit. 
        corrected_flux = flux/linear_regression_model   

    else:
        #print('No regressors used.')
        corrected_flux = flux
    
    return corrected_flux

def simple_lightcurve(short_name, phot_type='aper'):
    """Makes a 'simple' lightcurve with each reference star weighted equally when creating the artificial comparison lightcurve. Makes a light curve csv for every photometry file of type phot_type.

    :param short_name: short name of the target
    :type target: str
    :param phot_type: type of photometry, 'aper' or 'psf', defaults to 'aper'
    :type phot_type: str, optional
    """
    print('\nRunning simple_lightcurve().\n')

    pines_path = pines_dir_check()

    outlier_tolerance = 0.2 #If a reference > outlier_tolerance of its values above sigma clipping threshold, mark it as bad. 
    #centroided_sources.columns = centroided_sources.keys().str.strip()

    centroid_path = pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv')
    if not os.path.exists(centroid_path):
        raise RuntimeError('Centroid csv does not exist for {}, you need to make one using centroider.'.format(short_name))
    centroided_sources = pines_log_reader(centroid_path)
    sources = get_source_names(centroided_sources)

    #Get list of photometry files for this target. 
    photometry_path = pines_path/('Objects/'+short_name+'/'+phot_type+'_phot/')
    analysis_path = pines_path/('Objects/'+short_name+'/analysis/simple_lc_analysis/')
    if not os.path.exists(analysis_path):
        os.mkdir(analysis_path)
    photometry_files = natsorted([x for x in photometry_path.glob('*.csv')])

    num_refs = len(sources) - 1
    
    #Loop over all photometry files in the aper_phot directory. 
    for i in range(len(photometry_files)):
        #Load in the photometry data. 
        if phot_type =='aper':
            aperture_radius = float(str(photometry_files[i]).split('_')[-3])
        
        phot_data = pines_log_reader(photometry_files[i])

        #Remove entries that have NaN's for flux values.
        for j in range(len(sources)):
            name = sources[j]
            phot_data[name+' Flux'] = phot_data[name+' Flux'].astype(float)
            phot_data[name+' Flux Error'] = phot_data[name+' Flux Error'].astype(float)

        #Get target interpolation warnings. 
        targ_interp_flags = np.array(phot_data[short_name+' Interpolation Flag'])

        #Get times of exposures. 
        times = np.array(phot_data['Time BJD TDB'])
        
        #Convert to datetimes for plotting purposes.
        dts = np.array([julian.from_jd(times[i], fmt='jd') for i in range(len(times))])
        
        #Get the target's flux and background
        targ_flux = np.array(phot_data[short_name+' Flux'])
        targ_flux_err = np.array(phot_data[short_name+' Flux Error'])

        #Get the reference stars' fluxes and backgrounds. 
        ref_flux = np.zeros((num_refs, len(phot_data)))
        ref_flux_err = np.zeros((num_refs, len(phot_data)))
        for j in range(0, num_refs):
            ref_flux[j,:] = phot_data['Reference '+str(j+1)+' Flux']
            ref_flux_err[j,:] = phot_data['Reference '+str(j+1)+' Flux Error']
            #Discard variable stars. 
            #values, clow, chigh = sigmaclip(ref_flux[j], low=2.5, high=2.5)
            # if (len(phot_data) - len(values)) > (int(outlier_tolerance * len(phot_data))):
            #     print('Have to add flagging bad refs.')

        closest_ref = np.where(abs(np.nanmean(ref_flux, axis=1)-np.nanmean(targ_flux)) == min(abs(np.nanmean(ref_flux, axis=1)-np.nanmean(targ_flux))))[0][0]

        #Split data up into individual nights. 
        night_inds = night_splitter(times)
        num_nights = len(night_inds)
        
        #if plot_mode == 'combined':
        #    fig, axis = plt.subplots(nrows=1, ncols=num_nights, figsize=(16, 5))
        
        #colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

        #Get the time range of each night. Set each plot panel's xrange according to the night with the longest time. Makes seeing potential variability signals easier. 
        night_lengths = np.zeros(num_nights)
        for j in range(num_nights):
            inds = night_inds[j]
            night_lengths[j] = times[inds][-1] - times[inds][0]
        longest_night = max(night_lengths)
        longest_night_hours = np.ceil(longest_night*24)
        global line_list, filename_list
        line_list = []
        filename_list = []
        flux_save = []
        for j in range(num_nights):
            j += 1

            inds = night_inds[j-1]
            filename_list.append(np.array([phot_data['Filename'][z] for z in inds]))
            alc = np.zeros(len(inds))

            #Normalize reference lightcurves
            #TODO: Each night should be normalized separately. 
            for k in range(num_refs):
                ref_flux[k][inds] = ref_flux[k][inds] / np.nanmedian(ref_flux[k][inds])

            for k in range(len(inds)):
                #Do a sigma clip on normalized references to avoid biasing median. 
                ###values, clow, chigh = sigmaclip(ref_flux[:,inds[k]][~np.isnan(ref_flux[:,inds[k]])], low=1.5, high=1.5) 
                ###alc[k] = np.median(values)
                avg, med, std = sigma_clipped_stats(ref_flux[:,inds[k]][~np.isnan(ref_flux[:,inds[k]])], sigma=1.5)
                alc[k] = med
            
            #Correct the target lightcurve using the alc. 
            alc = alc / np.nanmedian(alc)
            targ_flux_norm = targ_flux[inds] / np.nanmedian(targ_flux[inds])
            targ_corr = targ_flux_norm / alc
            targ_corr = targ_corr / np.nanmedian(targ_corr)
            flux_save.extend(targ_corr)

            #Correct the example reference lightcurve using the alc. 
            ref_corr_norm = ref_flux[closest_ref][inds] / np.nanmedian(ref_flux[closest_ref][inds])
            ref_corr = ref_corr_norm / alc
            ref_corr = ref_corr / np.nanmedian(ref_corr)

            # #Do sigma clipping on the corrected lightcurve to get rid of outliers (from clouds, bad target centroid, cosmic rays, etc.)
            # ###vals, lo, hi = sigmaclip(targ_corr, low=2.5, high=2.5)
            # avg, med, std = sigma_clipped_stats(targ_corr, sigma=3)      
            # bad_vals = np.where((targ_corr > med + 5*std) | (targ_corr < med - 5*std))[0]
            # good_vals = np.where((targ_corr < med + 5*std) & (targ_corr > med - 5*std))[0]      

        #Output the simple lc data to a csv. 
        time_save = times        
        flux_err_save = np.zeros(len(flux_save)) + np.nanstd(targ_corr)
        output_dict = {'Time':time_save, 'Flux':flux_save, 'Flux Error':flux_err_save}
        output_df = pd.DataFrame(data=output_dict)
        if phot_type == 'aper':
            output_filename = analysis_path/(short_name+'_simple_lc_aper_phot_'+str(np.round(aperture_radius,1))+'_pix.csv')
            print('\nSaving to {}.\n'.format(output_filename))
            output_df.to_csv(output_filename)
        elif phot_type == 'psf':
            print("ERROR: Need to create flux output for PSF photometry.")     
               
    return

def weighted_lightcurve(short_name, phot_type='aper', convergence_threshold=1e-9, mode='night', n_sig_refs=5, sigma_clip_threshold=4., max_iterations=1000, use_pwv=False, red_stars_only=False, blue_stars_only=False, high_correlation_refs=False, force_output_path=''):
    """Creates a corrected light curve for the target using a weighted mean of reference star fluxes. 

    :param short_name: short name of the target
    :type short_name: str
    :param phot_type: photometry type, 'aper' or 'psf', defaults to 'aper'
    :type phot_type: str, optional
    :param convergence_threshold: threshold in reference star weighting loop to consider weights converged, defaults to 1e-9
    :type convergence_threshold: float, optional
    :param mode: 'night' or 'global' normalization mode, defaults to 'night'
    :type mode: str, optional
    :param n_sig_refs: n brightest reference stars whose average standard deviation is evaluated to choose the best aperture, defaults to 5
    :type n_sig_refs: int, optional
    :param sigma_clip_threshold: threshold for sigma clipping of target light curve, defaults to 4.0
    :type sigma_clip_threshold: float, optional
    :param max_iterations: maximum number of iterations allowed in reference star weighting loop, defaults to 1000
    :type max_iterations: int, optional
    :param use_pwv: whether or not to use precipitable water vapor data in the regression, defaults to False
    :type use_pwv: bool, optional
    :param red_stars_only: whether or not to use only the reddest stars in the field as references, defaults to False
    :type red_stars_only: bool, optional
    :param blue_stars_only: whether or not to use only the bluest stars in the field as references, defaults to False
    :type blue_stars_only: bool, optional
    :param high_correlation_refs: whether or not to use only the references that are most highly correlated with the raw target flux, defaults to False
    :type high_correlation_refs: bool, optional
    :param force_output_path: user-chosen path to use if you do not want to use the default ~/Documents/PINES_analysis_toolkit directory for analysis, defaults to ''
    :type force_output_path: str, optional
    """
    print('\nRunning weighted_lightcurve() in {} normalization mode.'.format(mode))

    if (phot_type != 'aper') and (phot_type != 'psf'):
        raise ValueError("phot_type must be either 'aper' or 'psf'!")
    
    def optimal_aperture_output():
        '''
            PURPOSE: writes the aperture type that minimizes the eference star scatter metric to a file 'optimal_aperture.txt' in the object's analysis directory.
        '''
        best_ap_output_path = analysis_path/('optimal_aperture.txt')
        if not os.path.exists(best_ap_output_path):
            with open(best_ap_output_path, 'w') as f:
                f.write('Nightly-normalized:  \n')
                f.write('Globally-normalized: ')

        #Write the best aperture to a text file.
        if mode == 'night':
            print('\nBest aperture for nightly-normalized photometry: {}.'.format(best_ap))
            with open(best_ap_output_path, 'r') as f:
                lines = f.readlines()
            
            lines[0] = 'Nightly-normalized:  '+best_ap+'\n'
            with open(best_ap_output_path, 'w') as f:
                for line in lines:
                    f.write(line)

        elif mode == 'global':
            print('\nBest aperture for globally-normalized photometry: {}.'.format(best_ap))
            with open(best_ap_output_path, 'r') as f:
                lines = f.readlines()
                    
            lines[1] = 'Globally-normalized: '+best_ap
            with open(best_ap_output_path, 'w') as f:
                for line in lines:
                    f.write(line)

    def optimizing_metric_caclulator(all_nights_binned_corr_ref_flux, n_sig_refs):
        all_ref_bin_avg_std = [[] for x in range(len(all_nights_binned_corr_ref_flux[0][0]))]
        for kk in range(num_nights):
            for jj in range(num_refs):
                all_ref_bin_avg_std[jj].extend([np.nanmean(np.nanstd(all_nights_binned_corr_ref_flux[kk][:,jj]))])

        #Use the average standard deviation of the binned data for the brightest n_sig_refs reference stars to determine the optimal aperture radius, fixed or variable. 
        #We use multiple reference stars (default = 5) to calculate this statistic, in case of any variability in the reference stars. 
        #We want to *minimize* binned_std to determine the best radius aperture. 
        return np.mean(all_ref_bin_avg_std[0:n_sig_refs])

    def output_filename_generator(analysis_path, ap_rad, mode):
        if phot_type == 'aper':
            if not os.path.exists(analysis_path/('aper_phot_analysis/'+ap_rad)):
                os.mkdir(analysis_path/('aper_phot_analysis/'+ap_rad))
            if 'fixed' in ap_rad:
                rad = ap_rad.split('_')[0]
                if mode == 'night':
                    output_filename = analysis_path/('aper_phot_analysis/'+ap_rad+'/'+short_name.replace(' ','')+'_fixed_aper_phot_r='+rad+'_nightly_weighted_lc.csv')
                elif mode == 'global':
                    output_filename = analysis_path/('aper_phot_analysis/'+ap_rad+'/'+short_name.replace(' ','')+'_fixed_aper_phot_r='+rad+'_global_weighted_lc.csv')

            elif 'variable' in ap_rad:
                factor = ap_rad.split('_')[0]
                if mode == 'night':
                    output_filename = analysis_path/('aper_phot_analysis/'+ap_rad+'/'+short_name.replace(' ','')+'_variable_aper_phot_f='+factor+'_nightly_weighted_lc.csv')
                elif mode == 'global':
                    output_filename = analysis_path/('aper_phot_analysis/'+ap_rad+'/'+short_name.replace(' ','')+'_variable_aper_phot_f='+factor+'_global_weighted_lc.csv')
        else:
            raise RuntimeError('Need to implement PSF photometry!')
        return output_filename

    np.seterr(divide='ignore') #Ignore divide-by-zero errors. 

    #Set up paths and get photometry files. 
    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()

    phot_path = pines_path/('Objects/'+short_name+'/'+phot_type+'_phot/')
    phot_files = natsorted(glob(str(phot_path/'*.csv')))
    analysis_path = pines_path/('Objects/'+short_name+'/analysis/')
    
    #Get source names for this observation. 
    source_detection_file = pines_path/('Objects/'+short_name+'/sources/target_and_references_source_detection.csv')
    source_df = pd.read_csv(source_detection_file)

    if red_stars_only and blue_stars_only:
        raise ValueError('red_stars_only and blue_stars_only cannot both be True!')

    #If red_stars_only, limit the source_df to the red stars.
    if red_stars_only:
        print('Using red reference stars only!')
        target_entry = source_df.loc[[0]]

        #Grab only the reddest sources from the sources DataFrame.
        source_df = source_df[(source_df['M_G'] < 14) & (source_df['M_G'] > 7.5)]

        #Create the new source_df with the target in the first row. 
        source_df = pd.concat([target_entry, source_df])

        #Reset the indices of the DataFrame.s
        source_df.reset_index(inplace=True, drop=True)

    if blue_stars_only:
        print('Using blue reference stars only!')
        target_entry = source_df.loc[[0]]

        #Grab only the reddest sources from the sources DataFrame.
        source_df = source_df[(source_df['M_G'] < 7.5)]

        #Create the new source_df with the target in the first row. 
        source_df = pd.concat([target_entry, source_df])

        #Reset the indices of the DataFrame.s
        source_df.reset_index(inplace=True, drop=True)
    
    source_names = np.array(source_df['Name'])
    ref_names = source_names[1:]
    num_refs = len(ref_names)

    #Get centroids for regression. 
    centroid_path = pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv')
    centroid_df = pines_log_reader(centroid_path)
    full_centroid_x = np.array(centroid_df[source_names[0]+' Image X'], dtype='float')
    full_centroid_y = np.array(centroid_df[source_names[0]+' Image Y'], dtype='float')

    if use_pwv:
        #Get PWV for regression. 
        pwv_path = pines_path/('Objects/'+short_name+'/pwv/'+short_name+'_fyodor_pwv.csv')
        pwv_df = pines_log_reader(pwv_path)
        full_pwv = np.array(pwv_df['PWV'])

    best_metric = 9999 #Initialize the metric that we'll use to choose the best lightcurve.

    #Loop over all photometry files.
    for i in range(len(phot_files)):
        phot_file = Path(phot_files[i])
        if phot_type == 'aper':
            if 'variable' in phot_file.name:
                ap_rad = phot_file.name.split('aper_phot_')[1].split('_')[0]+'_variable'
                print('\nRestoring variable aperture photometry output, seeing multiplicative factor = {}.'.format(phot_file.name.split('_')[4]))

            elif 'fixed' in phot_file.name:
                ap_rad = phot_file.name.split('aper_phot_')[1].split('_')[0]+'_fixed'
                print('\nRestoring fixed aperture photometry output, pixel radius = {}.'.format(ap_rad.split('_')[0]))
        else:
            ap_rad == ''

        df = pd.read_csv(phot_file)
        df.columns = df.keys().str.strip()
        
        full_times = np.array(df['Time BJD TDB'])        
        full_file_list = np.array(df['Filename'])
        full_airmass = np.array(df['Airmass']) #Get airmass for regression.
        full_background = np.array(df[source_names[0]+' Background'], dtype='float')
        sigma_clip_flags = np.zeros(len(full_times), dtype='int') #Track any frames identified as bad in sigma clipping. 

        #Split data up into individual nights. 
        if mode == 'night':
            night_inds = night_splitter(full_times)
        #Or treat all observations as a single "night".
        elif mode == 'global':
            night_inds = [list(np.arange(len(full_times)))]
        else:
            raise ValueError("mode must be either 'night' or 'global'!")

        num_nights = len(night_inds)

        #Declare some lists to hold data for all nights after processing is complete. 
        all_nights_file_list             = []
        all_nights_times                 = []
        all_nights_binned_times          = []

        all_nights_raw_targ_flux         = []
        all_nights_raw_targ_err          = []
        all_nights_norm_targ_flux        = []
        all_nights_norm_targ_err         = []

        all_nights_raw_ref_flux          = []
        all_nights_raw_ref_err           = []
        all_nights_norm_ref_flux         = []
        all_nights_norm_ref_err          = []

        all_nights_alc_flux              = []
        all_nights_alc_err               = []

        all_nights_corr_targ_flux        = []
        all_nights_corr_targ_err         = []
        all_nights_binned_corr_targ_flux = []
        all_nights_binned_corr_targ_err  = []

        all_nights_corr_ref_flux         = []
        all_nights_corr_ref_err          = []
        all_nights_binned_corr_ref_flux  = []
        all_nights_binned_corr_ref_err   = []

        norm_night_weights = []

        #Loop over each night in the dataset.
        for j in range(num_nights):
            num_frames = len(night_inds[j])
            inds = np.array(night_inds[j])
            num_frames = len(inds)
            airmass = full_airmass[inds]
            centroid_x = full_centroid_x[inds]
            centroid_y = full_centroid_y[inds]
            background = full_background[inds]
            if use_pwv:
                pwv = full_pwv[inds]     
                        
            times = np.array(full_times[inds])
            all_nights_times.append(times)
            date = julian.from_jd(times[0])
            date_str = 'UT '+date.strftime('%b %d, %Y')

            #Grab the target flux and error. 
            targ_flux = np.array(df[source_names[0]+' Flux'][inds], dtype='float')
            targ_flux_err = np.array(df[source_names[0]+' Flux Error'][inds], dtype='float')
            
            #Normalize the target flux and error. 
            normalization = np.nansum(targ_flux/targ_flux_err**2) / np.nansum(1/targ_flux_err**2)
            targ_flux_norm = targ_flux/normalization 

            #Sigmaclip the normalized target flux, in order to toss any suspect measurements. 
            vals, lower, upper = sigmaclip(targ_flux_norm, low=sigma_clip_threshold, high=sigma_clip_threshold)
            good_frames = np.where((targ_flux_norm > lower) & (targ_flux_norm < upper))[0]
            bad_frames = np.where((targ_flux_norm < lower) | (targ_flux_norm > upper))[0]
                
            #inds = inds[good_frames]
            bad_inds = inds[bad_frames]
            sigma_clip_flags[bad_inds] = 1 #Flag the frame as bad (1) for the output. 

            #Regrab raw flux using *good* frames. 
            targ_flux[bad_frames] = np.nan
            targ_flux_err[bad_frames] = np.nan

            #targ_flux = np.array(df[short_name+' Flux'][inds])
            #targ_flux_err = np.array(df[short_name+' Flux Error'][inds])
            all_nights_raw_targ_flux.append(targ_flux)
            all_nights_raw_targ_err.append(targ_flux_err)
          
            #Recompute the normalization using only the good frames. 
            ##good_frame_normalization = np.nansum(targ_flux/targ_flux_err**2) / np.nansum(1/targ_flux_err**2)
            good_frame_normalization = np.nanmean(targ_flux)
            targ_flux_norm = targ_flux / good_frame_normalization
            targ_flux_err_norm = targ_flux_err / good_frame_normalization
            all_nights_norm_targ_flux.append(targ_flux_norm)
            all_nights_norm_targ_err.append(targ_flux_err_norm)

            #make numpy arrays for all the stuff we need
            raw_flux =       np.zeros((num_frames,num_refs)) #units of photons
            raw_err =        np.zeros((num_frames,num_refs)) #units of photons
            norm_flux =      np.zeros((num_frames,num_refs)) #raw flux normalized by a single value to be near 1.0, unitless
            norm_err =       np.zeros((num_frames,num_refs)) #norm flux err, unitless
            

            try:    
                block_inds = block_splitter(times) #Get indices of each block for binning/plotting purposes.
            except:
                pdb.set_trace()
                
            files = np.array(full_file_list[inds])
            all_nights_file_list.append(files)

            #Write out bad frames, if any. 
            if len(bad_frames) > 0:
                print('Sigma clipping identified {} bad frames:'.format(len(bad_frames)))
                for f in range(len(bad_frames)):
                    print('     '+full_file_list[bad_inds[f]])

            #fill raw flux and err flux arrays, and normalize
            for k in np.arange(num_refs):
                r = ref_names[k]
                raw_flux[:,k] = np.array(df[r+' Flux'][inds], dtype='float')
                raw_err[:,k] = np.array(df[r+ ' Flux Error'][inds], dtype='float')
                raw_flux[bad_frames, k] = np.nan #Overwrite any sigma-clipped frames with NaNs. 
                raw_err[bad_frames, k] = np.nan
                
                #Normalize to weighted mean of flux values so they are near 1.
                #Weights are 1/raw_err**2.
                normalization = np.nansum(raw_flux[:,k]/raw_err[:,k]**2) / np.nansum(1/raw_err[:,k]**2) #ORIGINAL VERSION!

                norm_flux[:,k] = raw_flux[:,k] / normalization
                norm_err[:,k] = raw_err[:,k] / normalization   

            all_nights_raw_ref_flux.append(raw_flux)
            all_nights_raw_ref_err.append(raw_err)
            all_nights_norm_ref_flux.append(norm_flux)
            all_nights_norm_ref_err.append(norm_err)

            if high_correlation_refs:
                selection_type = 'worst'
                n_to_take = 4

                #Measure the Pearson correlation coefficient of each reference normalized flux with the target's normalized flux.
                correlations = np.zeros(num_refs)
                for k in range(num_refs):
                    non_nan_inds = np.where(~np.isnan(norm_flux[:,k]) | ~np.isnan(targ_flux_norm))[0]
                    corr, sig = pearsonr(targ_flux_norm[non_nan_inds], norm_flux[:,k][non_nan_inds])
                    correlations[k] = corr

                #Sort the references by those with the highest correlation with the target, and take the best handful of them. 
                sorted_corr_inds = np.argsort(correlations)[::-1]

                #Can choose the n_to_take most correlated references...
                if selection_type == 'best':
                    best_ref_inds = sorted_corr_inds[0:n_to_take]
                #Or the least correlated, for testing purposes...
                elif selection_type == 'worst':
                    best_ref_inds = sorted_corr_inds[n_to_take:]

                #Change the norm_flux array to only contain the flux from the most highly correlated references.
                norm_flux = norm_flux[:,best_ref_inds]
                norm_err  = norm_err[:,best_ref_inds]
                num_refs = len(norm_flux[0,:])
                ref_names = ref_names[best_ref_inds]

            #Define some more arrays. 
            alc_flux =       np.zeros((num_frames,num_refs)) #special ALC to remove atm veriation, unitless, near 1.0
            alc_err =        np.zeros((num_frames,num_refs)) #err in special ALC
            corr_flux =      np.zeros((num_frames,num_refs)) #norm flux corrected to remove atm veriation, unitless, near 1.0
            corr_err =       np.zeros((num_frames,num_refs)) #corr flux err
            regressed_corr_flux = np.zeros((num_frames,num_refs)) #corrected flux that has been run through the linear regression procedure

            #Now calculate the "special" ALC for each ref star, use that to correct the ref star's flux, and run the corrected flux through the linear regression model.
            #Regressed corrected flux is used to measure the standard deviation of the lightcurve, which is used to weight it in the target ALC. 
            old_stddev = np.zeros(num_refs)
            new_stddev = np.zeros(num_refs)
            keep_going = True
            iteration = 0 
            while keep_going:
                iteration += 1

                for k in np.arange(num_refs):
                    #indices w/o the k'th ref star
                    indices = np.where(np.arange(num_refs) != k)
                    
                    #make local arrays to store non_k norm fluxes
                    norm_flux_without_k = np.resize(norm_flux[:,indices],(num_frames,num_refs-1))
                    norm_err_without_k = np.resize(norm_err[:,indices],(num_frames,num_refs-1))

                    #ALC = weighted mean of norm flux, weights = 1/norm_err**2
                    alc_flux[:,k] = np.sum(norm_flux_without_k/norm_err_without_k**2,axis=1) / np.sum(1/norm_err_without_k**2,axis=1)
                    alc_err[:,k] = np.sqrt( 1 / np.sum(1/norm_err_without_k**2,axis=1))
                          
                    #Renormalize the ALC and normalized reference flux so that their mean is EXACTLY 1. 
                    alc_renorm = np.nanmean(alc_flux[:,k])
                    alc_flux[:,k] = alc_flux[:,k] / alc_renorm
                    norm_renorm = np.nanmean(norm_flux[:,k])
                    norm_flux[:,k] = norm_flux[:,k] / norm_renorm

                    #calculate corr flux and corr err
                    corr_flux[:,k] = norm_flux[:,k] / alc_flux[:,k]
                    corr_err[:,k] = corr_flux[:,k] * np.sqrt((norm_err[:,k]/norm_flux[:,k])**2 + (alc_err[:,k]/alc_flux[:,k])**2)

                    #Perform a regression of the corrected flux against various parameters. 
                    if use_pwv:
                        regression_dict = {'airmass':airmass, 'centroid_x':centroid_x, 'centroid_y':centroid_y, 'pwv':pwv}
                    else:
                        regression_dict = {'airmass':airmass, 'centroid_x':centroid_x, 'centroid_y':centroid_y}

                    regressed_corr_flux[:,k] = regression(corr_flux[:,k], regression_dict)

                    #Calculate stddevs from the regressed corrected flux.
                    old_stddev[k] = new_stddev[k]
                    new_stddev[k] = np.nanstd(regressed_corr_flux[:,k])

                    #update normalized errors
                    norm_err[:,k] = new_stddev[k]

                fractional_changes = np.abs(new_stddev - old_stddev)/old_stddev
                
                #Once fractional changes stabilize, break out of the loop. 
                if np.sum(fractional_changes) < convergence_threshold:
                    keep_going = False

                if iteration > max_iterations:
                    print('Iterations in reference star weighting loop exceeded max_iterations.')
                    break
                
            
            all_nights_corr_ref_flux.append(regressed_corr_flux)
            all_nights_corr_ref_err.append(corr_err)

            alc_final_flux = np.sum(norm_flux/norm_err**2,axis=1) / np.sum(1/norm_err**2,axis=1)
            alc_final_err = np.sqrt( 1 / np.sum(1/norm_err**2,axis=1))    

            norm_night_weights.append((1/norm_err[0,:])**2/np.sum((1/norm_err[0,:])**2)) #Save the NORMALIZED weight for each comparison star on this night.
            all_nights_alc_flux.append(alc_final_flux)
            all_nights_alc_err.append(alc_final_err)

            #Correct the target flux with the master ALC and renormalize to 1. 
            targ_flux_corr = targ_flux_norm/alc_final_flux
            targ_flux_corr_err = np.sqrt( (targ_flux_err_norm/alc_final_flux)**2 + (targ_flux_norm*alc_final_err/(alc_final_flux**2))**2)
            
            if use_pwv:
                regression_dict = {'airmass':airmass, 'centroid_x':centroid_x, 'centroid_y':centroid_y,'pwv':pwv}
            else: 
                regression_dict = {'airmass':airmass, 'centroid_x':centroid_x, 'centroid_y':centroid_y}
            
            regressed_targ_flux_corr = regression(targ_flux_corr, regression_dict, verbose=True)

            
            #regressed_targ_flux_corr = leave_one_out_regression(times, targ_flux_corr, targ_flux_corr_err, regression_dict, verbose=False)
            #pdb.set_trace()

            all_nights_corr_targ_flux.append(regressed_targ_flux_corr)
            all_nights_corr_targ_err.append(targ_flux_corr_err)

            binned_times = np.zeros(len(block_inds))
            binned_flux  = np.zeros(len(block_inds))
            binned_errs  = np.zeros(len(block_inds))

            binned_ref_fluxes = np.zeros((len(block_inds), num_refs))
            binned_ref_flux_errs = np.zeros((len(block_inds), num_refs))

            for k in range(len(block_inds)):
                binned_times[k] = np.nanmean(times[block_inds[k]])
                binned_flux[k]  = np.nanmean(regressed_targ_flux_corr[block_inds[k]])
                #binned_errs[k] = np.nanmean(targ_flux_corr_err[block_inds[k]]) / np.sqrt(len(block_inds[k])) #Uses calculated errors. 
                binned_errs[k]  = np.nanstd(targ_flux_corr[block_inds[k]])/np.sqrt(len(block_inds[k])) #Uses standard deviation of the block.
                
                for l in range(num_refs):
                    binned_ref_fluxes[k,l] = np.nanmean(regressed_corr_flux[:,l][block_inds[k]])
                   #binned_ref_flux_errs[k,l] = np.nanmean(corr_err[:,l][block_inds[k]]) / np.sqrt(len(block_inds[k])) #Uses calculated errors.
                    binned_ref_flux_errs[k,l] = np.std(corr_flux[:,l][block_inds[k]])/np.sqrt(len(block_inds[k])) #Uses standard deviation of the block.

            all_nights_binned_times.append(binned_times)
            all_nights_binned_corr_targ_flux.append(binned_flux)
            all_nights_binned_corr_targ_err.append(binned_errs)
            all_nights_binned_corr_ref_flux.append(binned_ref_fluxes)
            all_nights_binned_corr_ref_err.append(binned_ref_flux_errs)

        #Calculate our metric for choosing our best lightcurve. 
        optimizing_metric = optimizing_metric_caclulator(all_nights_binned_corr_ref_flux, n_sig_refs)

        print('Average stddev of binned corrected reference flux for brightest {} references: {:1.5f}'.format(n_sig_refs, optimizing_metric))
        if optimizing_metric < best_metric:
            print('New best lightcurve found!')
            best_metric = optimizing_metric
            best_ap = ap_rad

        #Write out the best aperture 
        optimal_aperture_output()
        
        #Output the weighted lc data to a csv. 
        #Have to unpack things that were saved into night arrays. 
        
        file_list_save = full_file_list
        time_save = full_times
        sigma_clip_flag_save = sigma_clip_flags

        targ_flux_corr_save = []
        targ_flux_corr_err_save = []
        targ_flux_norm_save = []
        targ_err_norm_save = []
        alc_flux_save = []
        alc_err_save = []
        ref_flux_corr_save = [[] for x in range(num_refs)]
        ref_flux_corr_err_save = [[] for x in range(num_refs)]
        night_weight_arrays = [[] for x in range(num_refs)]
        for j in range(num_nights):
            targ_flux_corr_save.extend(all_nights_corr_targ_flux[j])
            targ_flux_corr_err_save.extend(all_nights_corr_targ_err[j])
            targ_flux_norm_save.extend(all_nights_norm_targ_flux[j])
            targ_err_norm_save.extend(all_nights_norm_targ_err[j])
            alc_flux_save.extend(all_nights_alc_flux[j])
            alc_err_save.extend(all_nights_alc_err[j])
            for k in range(num_refs):
                ref_flux_corr_save[k].extend(all_nights_corr_ref_flux[j][:,k])
                ref_flux_corr_err_save[k].extend(all_nights_corr_ref_err[j][:,k])
                night_weight_arrays[k].extend(np.zeros(len(night_inds[j]))+norm_night_weights[j][k])
        
        #Create the output dictionary
        output_dict = {'Filename':file_list_save, 'Time BJD TDB':time_save, 'Sigma Clip Flag':sigma_clip_flag_save, short_name+' Normalized Flux':targ_flux_norm_save, short_name+' Normalized Flux Error':targ_err_norm_save, short_name+' Corrected Flux':targ_flux_corr_save, short_name+' Corrected Flux Error':targ_flux_corr_err_save, 'ALC Flux':alc_flux_save, 'ALC Flux Error':alc_err_save}
        for j in range(num_refs):
            ref = ref_names[j]
            output_dict[ref+' Corrected Flux'] = ref_flux_corr_save[j]
            output_dict[ref+' Corrected Flux Error'] = ref_flux_corr_err_save[j]
        
        #Convert to a dataframe
        output_df = pd.DataFrame(data=output_dict)
        output_filename = output_filename_generator(analysis_path, ap_rad, mode)
        print('Saving weighted_lightcurve output to {}.'.format(output_filename))

        #Write out to a csv file.
        with open(output_filename, 'w') as f:
            for j in range(len(output_df)):
                #Write the header line. 
                if j == 0:
                    f.write('{:<21s}, '.format('Filename'))
                    f.write('{:<15s}, '.format('Time BJD TDB'))
                    f.write('{:<15s}, '.format('Sigma Clip Flag'))
                    f.write('{:<31s}, {:<31s}, {:<30s}, {:<30s}, {:<9s}, {:<15s}, '.format(short_name+' Normalized Flux',short_name+' Normalized Flux Error',short_name+' Corrected Flux', short_name+' Corrected Flux Error', 'ALC Flux','ALC Flux Error'))

                    for i in range(len(ref_names)):
                        n = ref_names[i]
                        if i != num_refs - 1:
                            f.write('{:<22s}, {:<30s}, {:<36s}, '.format(n+' ALC Weight', n+' Corrected Flux', n+' Corrected Flux Error'))
                        else:
                            f.write('{:<22s}, {:<30s}, {:<36s}\n'.format(n+' ALC Weight', n+' Corrected Flux', n+' Corrected Flux Error'))
                
                #Write in the data lines.
                f.write('{:<21s}, '.format(output_df['Filename'][j]))
                f.write('{:<15.7f}, '.format(output_df['Time BJD TDB'][j]))
                f.write('{:<15d}, '.format(output_df['Sigma Clip Flag'][j]))
                f.write('{:<31.6f}, {:<31.6f}, {:<30.6f}, {:<36.6f}, {:<9.6f}, {:<15.6f}, '.format(output_df[short_name+' Normalized Flux'][j],output_df[short_name+' Normalized Flux Error'][j],output_df[short_name+' Corrected Flux'][j], output_df[short_name+' Corrected Flux Error'][j], output_df['ALC Flux'][j], output_df['ALC Flux Error'][j]))
                
                for i in range(num_refs):
                    n = ref_names[i]                  
                    if i != num_refs - 1:
                        format_string = '{:<22.3f}, {:<30.6f}, {:<36.6f}, '
                    else:
                        format_string = '{:<22.3f}, {:<30.6f}, {:<36.6f}\n'
                    f.write(format_string.format(night_weight_arrays[i][j],output_df[n+' Corrected Flux'][j], output_df[n+' Corrected Flux Error'][j]))
    print('\nBest aperture: {}.\n'.format(best_ap))
    time.sleep(2)
    return
