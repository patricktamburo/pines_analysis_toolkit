import matplotlib.pyplot as plt 
import numpy as np 
import julian 
from astropy.stats import sigma_clipped_stats
import os
from pines_analysis_toolkit.analysis.block_binner import block_binner
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader
from pines_analysis_toolkit.utils.get_source_names import get_source_names
from pines_analysis_toolkit.analysis.night_splitter import night_splitter
from pines_analysis_toolkit.analysis.block_splitter import block_splitter

plt.ioff() 
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
    