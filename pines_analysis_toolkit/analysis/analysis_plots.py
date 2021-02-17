import matplotlib.pyplot as plt 
import numpy as np 
import pdb
import julian 
import matplotlib.dates as mdates
from astropy.stats import sigma_clipped_stats
import os 
plt.ioff() 
def standard_x_range(times):
    '''
        Calculates a standard range of time to use for x range in the plots for this target. 
    '''
    time_durations = np.zeros(len(times))
    for i in range(len(times)):
        time_durations[i] = (times[i][-1] - times[i][0])*24 + 0.25 #Gets the span of the night's observations in **HOURS**. Adds on an additional .25 hours for padding.
        
    return np.ceil(np.max(time_durations))/24 #Returns maximum time span of the nights in **DAYS** 

def standard_y_range(y, multiple=3.5):
    '''
        Calculates a standard y range for this target using the multiple * the standard deviation of the unbinned data.
    '''
    stds = np.zeros(len(y))
    for i in range(len(y)):
        stds[i] = sigma_clipped_stats(y[i])[2]
    return multiple*np.max(stds)

def raw_flux_plot(times, raw_targ_flux, raw_targ_err, raw_ref_flux, raw_ref_err, short_name, analysis_path, phot_type, ap_rad): 
    '''Authors:
        Patrick Tamburo, Boston University, December 2020
	Purpose:
        Creates a standard plot of nightly raw flux for target and reference stars. 
	Inputs:
        times (list of numpy arrays): Times of exposures (as Julian dates) for each night of observation.
        raw_targ_flux (list of numpy arrays): Raw flux for the target on each night of observation.
        raw_targ_err (list of numpy arrays): Flux errors for the target on each night of observation. 
        raw_ref_flux (list of numpy arrays): Raw flux for the reference stars on each night of observation. (n nights x n_exp per night x n_refss)
        raw_ref_err (list of numpy arrays): Flux errors for the reference stars on each night of observation.
        short_name (str): The object's short name.
        analysis_path (pathlib Path object): The path to this object's analysis directory.
    Outputs:

	TODO:

    '''
    
    num_nights = len(times)
    num_refs = np.shape(raw_ref_flux[0])[1]

    cm = plt.get_cmap('viridis_r')
    markers = ['+', 'x', '*', 'X']

    standard_time_range = standard_x_range(times)

    fig, axis = plt.subplots(nrows=1, ncols=num_nights, figsize=(17,5), sharey=True)
    
    for i in range(num_nights):
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
            ax.errorbar(times[i], raw_ref_flux[i][:,j], raw_ref_err[i][:,j], marker=marker, linestyle='-', label='Ref. '+str(j+1), color=color)
        
        #Set labels. 
        if i == 0:
            ax.set_ylabel('Photons', fontsize=20)
        ax.tick_params(labelsize=16)
        ax.set_xlabel('Time (JD$_{UTC}$)', fontsize=20)

        #Plot the target. 
        ax.errorbar(times[i], raw_targ_flux[i], raw_targ_err[i], marker='o', linestyle='-', label='Target', mfc='none', mew=2, color='tab:orange')

        if i == num_nights - 1:
            ax.legend(bbox_to_anchor=(1.03, 1), fontsize=14)
        
        ax.set_xlim(np.mean(times[i])-standard_time_range/2, np.mean(times[i])+standard_time_range/2)

    plt.suptitle(short_name+' Nightly Raw Flux', fontsize=20)
    plt.subplots_adjust(left=0.07, wspace=0.05, top=0.92, bottom=0.17)

    print('Saving raw flux figure...')
    if phot_type == 'aper':
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_'+'r='+ap_rad+'_nightly_raw_flux.png'
        if not os.path.exists(analysis_path/('aper_phot_analysis/'+ap_rad)):
            os.mkdir(analysis_path/('aper_phot_analysis/'+ap_rad))
        output_path = analysis_path/('aper_phot_analysis/'+ap_rad+'/'+filename)
        plt.savefig(output_path, dpi=300)
    else:
        raise RuntimeError("HAVE TO UPDATE!!!")
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_nightly_raw_flux.png'


def global_raw_flux_plot(times, raw_targ_flux, raw_targ_err, raw_ref_flux, raw_ref_err, short_name, analysis_path, phot_type, ap_rad):
    '''Authors:
        Patrick Tamburo, Boston University, December 2020
	Purpose:
        Creates a standard plot of global raw flux for target and reference stars. 
	Inputs:
        times (list of numpy arrays): Times of exposures (as Julian dates) for each night of observation.
        raw_targ_flux (list of numpy arrays): Raw flux for the target on each night of observation.
        raw_targ_err (list of numpy arrays): Flux errors for the target on each night of observation. 
        raw_ref_flux (list of numpy arrays): Raw flux for the reference stars on each night of observation. (n nights x n_exp per night x n_refss)
        raw_ref_err (list of numpy arrays): Flux errors for the reference stars on each night of observation.
        short_name (str): The object's short name.
        analysis_path (pathlib Path object): The path to this object's analysis directory.
    Outputs:

	TODO:

    '''
    times = times[0]
    num_refs = np.shape(raw_ref_flux[0])[1]

    cm = plt.get_cmap('viridis_r')
    markers = ['+', 'x', '*', 'X']

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(17,5), sharey=True)
    ax.set_yscale('log')
    for j in range(num_refs):
        color = cm(1.*j/num_refs) 
        marker = markers[j % len(markers)]
        ax.errorbar(times, raw_ref_flux[0][:,j], raw_ref_err[0][:,j], marker=marker, linestyle='-', label='Ref. '+str(j+1), color=color)
    ax.set_ylabel('Photons', fontsize=20)
    ax.tick_params(labelsize=16)
    ax.set_xlabel('Time (JD$_{UTC}$)', fontsize=20)

    ax.errorbar(times, raw_targ_flux[0], raw_targ_err[0], marker='o', linestyle='-', label='Target', mfc='none', mew=2, color='tab:orange')

    ax.legend(bbox_to_anchor=(1.01, 1), fontsize=14)
        
    plt.suptitle(short_name+' Global Raw Flux', fontsize=20)
    plt.subplots_adjust(left=0.07, wspace=0.05, top=0.92, bottom=0.17)

    print('Saving raw flux figure...')
    if phot_type == 'aper':
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_'+'r='+ap_rad+'_global_raw_flux.png'
        if not os.path.exists(analysis_path/('aper_phot_analysis/'+ap_rad)):
            os.mkdir(analysis_path/('aper_phot_analysis/'+ap_rad))
        output_path = analysis_path/('aper_phot_analysis/'+ap_rad+'/'+filename)
        plt.savefig(output_path, dpi=300)
    else:
        raise RuntimeError("HAVE TO UPDATE!!!")
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_nightly_raw_flux.png'

def normalized_flux_plot(times, norm_targ_flux, norm_targ_err, norm_ref_flux, norm_ref_err, alc_flux, alc_err, short_name, analysis_path, phot_type, ap_rad):
    '''Authors:
            Patrick Tamburo, Boston University, December 2020
        Purpose:
            Creates a standard plot of nightly normalized flux for target and reference stars. 
        Inputs:
            times (list of numpy arrays): Times of exposures (as Julian dates) for each night of observation.
            norm_targ_flux (list of numpy arrays): Normalized flux for the target on each night of observation.
            norm_targ_err (list of numpy arrays): Normalized flux errors for the target on each night of observation. 
            norm_ref_flux (list of numpy arrays): Normalized flux for the reference stars on each night of observation. (n nights x n_exp per night x n_refss)
            norm_ref_err (list of numpy arrays): Normalized flux errors for the reference stars on each night of observation.
            alc_flux (list of numpy arrays): The artficial comparison lightcurve on each night of observation.
            short_name (str): The object's short name.
            analysis_path (pathlib Path object): The path to this object's analysis directory.
        Outputs:

        TODO:

        '''
    
    num_nights = len(times)
    num_refs = np.shape(norm_ref_flux[0])[1]

    #cm = plt.get_cmap('viridis_r')
    #markers = ['+', 'x', '*', 'X']

    standard_time_range = standard_x_range(times)
    standard_y = standard_y_range(norm_targ_flux, multiple=3.0)

    fig, axis = plt.subplots(nrows=1, ncols=num_nights, figsize=(17,5), sharey=True)
   
    for i in range(num_nights):
        #Set the axis behavior depending if there are 1 or more nights. 
        if num_nights == 1:
            ax = axis
        elif num_nights > 1:
            ax = axis[i]

        #Plot the reference stars. 
        #for j in range(num_refs):
            #color = cm(1.*j/num_refs) 
            #marker = markers[j % len(markers)]
            # if j == 0:
            #     ax.errorbar(times[i], norm_ref_flux[i][:,j], norm_ref_err[i][:,j], marker='.', linestyle='', label='Ref. flux', color='tab:blue', zorder=0, alpha=0.2)
            # else:
            #     ax.errorbar(times[i], norm_ref_flux[i][:,j], norm_ref_err[i][:,j], marker='.', linestyle='', color='tab:blue', zorder=0, alpha=0.2)
        
        if i == 0:
            ax.set_ylabel('Normalized Flux', fontsize=20)
        ax.tick_params(labelsize=16)
        ax.set_xlabel('Time (JD$_{UTC}$)', fontsize=20)

        #Plot the target. 
        ax.errorbar(times[i], norm_targ_flux[i], norm_targ_err[i], color='tab:orange', mfc='none', marker='.', ms=7, lw=1.5, linestyle='', label='Target', mew=2, zorder=1)
        
        #Plot the ALC. 
        ax.errorbar(times[i], alc_flux[i], alc_err[i], color='tab:blue', lw=2, marker='.', linestyle='', label='ALC', zorder=2, mfc='none', mew=2, ms=5)
        
        if i == num_nights - 1:
            ax.legend(bbox_to_anchor=(1.03, 1), fontsize=14)
        
        ax.set_xlim(np.mean(times[i])-standard_time_range/2, np.mean(times[i])+standard_time_range/2)
        ax.set_ylim(1-standard_y, 1+standard_y)
        ax.grid(alpha=0.2)

    plt.suptitle(short_name+' Nightly Normalized Flux', fontsize=20)
    plt.subplots_adjust(left=0.07, wspace=0.05, top=0.92, bottom=0.17)

    #Save the figure. 
    print('Saving normalized flux figure...')
    if phot_type == 'aper':
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_'+'r='+ap_rad+'_nightly_normalized_flux.png'
        if not os.path.exists(analysis_path/('aper_phot_analysis/'+ap_rad)):
            os.mkdir(analysis_path/('aper_phot_analysis/'+ap_rad))
        output_path = analysis_path/('aper_phot_analysis/'+ap_rad+'/'+filename)
        plt.savefig(output_path, dpi=300)
    else:
        raise RuntimeError("HAVE TO UPDATE!!!")
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_nightly_raw_flux.png'
        
def global_normalized_flux_plot(times, norm_targ_flux, norm_targ_err, norm_ref_flux, norm_ref_err, alc_flux, alc_err, short_name, analysis_path, phot_type, ap_rad):
    '''Authors:
            Patrick Tamburo, Boston University, December 2020
        Purpose:
            Creates a standard plot of global normalized flux for target and reference stars. 
        Inputs:
            times (list of numpy arrays): Times of exposures (as Julian dates) for each night of observation.
            norm_targ_flux (list of numpy arrays): Normalized flux for the target on each night of observation.
            norm_targ_err (list of numpy arrays): Normalized flux errors for the target on each night of observation. 
            norm_ref_flux (list of numpy arrays): Normalized flux for the reference stars on each night of observation. (n nights x n_exp per night x n_refss)
            norm_ref_err (list of numpy arrays): Normalized flux errors for the reference stars on each night of observation.
            alc_flux (list of numpy arrays): The artficial comparison lightcurve on each night of observation.
            short_name (str): The object's short name.
            analysis_path (pathlib Path object): The path to this object's analysis directory.
        Outputs:

        TODO:

        '''
    times = times[0]
    num_refs = np.shape(norm_ref_flux[0])[1]

    #cm = plt.get_cmap('viridis_r')
    #markers = ['+', 'x', '*', 'X']

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(17,5))
    
    #for j in range(num_refs):
    #    color = cm(1.*j/num_refs) 
    #    marker = markers[j % len(markers)]
    #    ax.errorbar(times, norm_ref_flux[0][:,j], norm_ref_err[0][:,j], marker=marker, linestyle='', label='Ref. '+str(j+1), color=color)
    
    ax.set_ylabel('Normalized Flux', fontsize=20)
    ax.tick_params(labelsize=16)
    ax.set_xlabel('Time (JD$_{UTC}$)', fontsize=20)

    #Plot the target.
    ax.errorbar(times, norm_targ_flux[0], norm_targ_err[0], color='tab:orange', mfc='none', marker='.', ms=7, lw=1.5, linestyle='', label='Target', mew=2, zorder=1)
    
    #Plot the ALC. 
    ax.errorbar(times, alc_flux[0], alc_err[0], color='tab:blue', lw=2, marker='.', linestyle='', label='ALC', zorder=2, mfc='none', mew=2, ms=5)

    ax.legend(bbox_to_anchor=(1.01, 1), fontsize=14)
        
    plt.suptitle(short_name+' Global Normalized Flux', fontsize=20)
    plt.subplots_adjust(left=0.07, wspace=0.055, top=0.92, bottom=0.17)

    #Save the figure. 
    print('Saving normalized flux figure...')
    if phot_type == 'aper':
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_'+'r='+ap_rad+'_global_normalized_flux.png'
        if not os.path.exists(analysis_path/('aper_phot_analysis/'+ap_rad)):
            os.mkdir(analysis_path/('aper_phot_analysis/'+ap_rad))
        output_path = analysis_path/('aper_phot_analysis/'+ap_rad+'/'+filename)
        plt.savefig(output_path, dpi=300)
    else:
        raise RuntimeError("HAVE TO UPDATE!!!")

def corr_target_plot(times, targ_flux_corr, binned_times, binned_flux, binned_errs, short_name, analysis_path, phot_type, ap_rad):
    '''Authors:
        Patrick Tamburo, Boston University, November, December 2020
	Purpose:
        Creates a standard plot of nightly corrected target flux. 
	Inputs:
        target (numpy array): times of exposures (as Julian dates).
        targ_flux_corr (numpy array): corrected target flux measurements.
        binned_times (numpy array): times of block midpoints.
        binned_flux (numpy array): average fluxes of blocks.
        binned_errs (numpy array): uncertainties on binned fluxes.
        analysis_path (pathlib Path): path to the object's analysis directory.
        time_mode (str, optional): either 'UT' or 'JD'. UT for calendar-like times, JD for Julian dates. 
    Outputs:

	TODO:

    '''
    short_name = str(analysis_path).split('/')[-2]

    num_nights = len(times)

    standard_time_range = standard_x_range(times)
    standard_y = standard_y_range(targ_flux_corr)

    #Plot the corrected target flux. 
    fig, axis = plt.subplots(nrows=1,ncols=num_nights,figsize=(17,5), sharey=True)
    plt.subplots_adjust(left=0.07, wspace=0.05, top=0.86, bottom=0.17, right=0.96)

    for i in range(num_nights):
        #Set the axis behavior depending if there are 1 or more nights. 
        if num_nights == 1:
            ax = axis
        elif num_nights > 1:
            ax = axis[i]

        ax.axhline(1, color='r', lw=2, zorder=0)
        ax.grid(alpha=0.4)
        #Plot the corrected target and corrected binned target flux. 
        ax.plot(times[i], targ_flux_corr[i], linestyle='', marker='o', zorder=1, color='darkgrey', ms=5)
        ax.errorbar(binned_times[i], binned_flux[i], binned_errs[i], linestyle='', marker='o', zorder=2, color='k', capsize=2, ms=7)
        
        #Plot the 5-sigma decrement line. 
        five_sig = 5*np.mean(binned_errs[i])
        ax.axhline(1-five_sig, color='k', lw=2, linestyle='--', zorder=0)
        ax.tick_params(labelsize=16)
        
        #Set labels.
        if i == 0:
            ax.set_ylabel('Normalized Flux', fontsize=20)
        ax.set_xlabel('Time (JD$_{UTC}$)', fontsize=20)
        ut_date = julian.from_jd(times[i][0])
        ut_date_str = 'UT '+ut_date.strftime('%b. %d %Y')
        ax.set_title(ut_date_str, fontsize=20)
        ax.set_xlim(np.mean(times[i])-standard_time_range/2, np.mean(times[i])+standard_time_range/2)
        ax.set_ylim(1-standard_y, 1+standard_y)
    #plt.suptitle(short_name+' Nightly Corrected Target Lightcurve', fontsize=20)
    plt.subplots_adjust(left=0.07, wspace=0.05, top=0.92, bottom=0.17)

    #Save the figure. 
    print('Saving target lightcurve...')
    if phot_type == 'aper':
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_'+'r='+ap_rad+'_nightly_target_flux.png'
        if not os.path.exists(analysis_path/('aper_phot_analysis/'+ap_rad)):
            os.mkdir(analysis_path/('aper_phot_analysis/'+ap_rad))
        output_path = analysis_path/('aper_phot_analysis/'+ap_rad+'/'+filename)
        plt.savefig(output_path, dpi=300)
    else:
        raise RuntimeError("HAVE TO UPDATE!!!")
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_nightly_raw_flux.png'

def global_corr_target_plot(times, targ_flux_corr, binned_times, binned_flux, binned_errs, short_name, analysis_path, phot_type, ap_rad):
    '''Authors:
        Patrick Tamburo, Boston University, November, December 2020
	Purpose:
        Creates a standard plot of global corrected target flux. 
	Inputs:
        target (numpy array): times of exposures (as Julian dates).
        targ_flux_corr (numpy array): corrected target flux measurements.
        binned_times (numpy array): times of block midpoints.
        binned_flux (numpy array): average fluxes of blocks.
        binned_errs (numpy array): uncertainties on binned fluxes.
        analysis_path (pathlib Path): path to the object's analysis directory.
        time_mode (str, optional): either 'UT' or 'JD'. UT for calendar-like times, JD for Julian dates. 
    Outputs:

	TODO:

    '''
    times = times[0]
    binned_times = binned_times[0]
    short_name = str(analysis_path).split('/')[-2]


    #Plot the corrected target flux. 
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(17,5))
    plt.subplots_adjust(left=0.07, wspace=0.05, top=0.86, bottom=0.17, right=0.96)

    ax.axhline(1, color='r', lw=2, zorder=0)
    ax.grid(alpha=0.4)
    ax.plot(times, targ_flux_corr[0], linestyle='', marker='o', zorder=1, color='darkgrey', ms=5, label='Unbinned data')
    ax.errorbar(binned_times, binned_flux[0], binned_errs[0], linestyle='', marker='o', zorder=2, color='k', capsize=2, ms=7, label='Binned data')
    five_sig = 5*np.mean(binned_errs[0])
    ax.axhline(1-five_sig, color='k', lw=2, linestyle='--', zorder=0, label='5-$\sigma$ threshold')
    ax.tick_params(labelsize=16)
    
    ax.set_ylabel('Normalized Flux', fontsize=20)

    ax.set_xlabel('Time (JD$_{UTC}$)', fontsize=20)
    plt.suptitle(short_name+' Global Corrected Target Lightcurve', fontsize=20)
    
    #Save the figure. 
    print('Saving target lightcurve...')
    if phot_type == 'aper':
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_'+'r='+ap_rad+'_global_target_flux.png'
        if not os.path.exists(analysis_path/('aper_phot_analysis/'+ap_rad)):
            os.mkdir(analysis_path/('aper_phot_analysis/'+ap_rad))
        output_path = analysis_path/('aper_phot_analysis/'+ap_rad+'/'+filename)
        plt.savefig(output_path, dpi=300)
    else:
        raise RuntimeError("HAVE TO UPDATE!!!")

def corr_all_sources_plot(times, targ_flux_corr, binned_times, binned_flux, binned_errs, corr_flux, binned_ref_fluxes, binned_ref_flux_errs, short_name, analysis_path, phot_type, ap_rad, num_refs, num_nights, norm_night_weights):
    cmap = plt.cm.viridis
    fig, axis = plt.subplots(num_refs+1, num_nights, figsize=(14,4*num_refs))
    plt.subplots_adjust(top=0.95, left=0.1, bottom=0.1, right=0.95, hspace=0.1)
    
    standard_time_range = standard_x_range(times)
    standard_y = standard_y_range(targ_flux_corr)

    if num_nights > 1:
        #First, plot the target's corrected flux for each night. 
        ax = axis[0,:]
        for i in range(num_nights):
            ax[i].axhline(1, color='k', alpha=0.7, lw=1, zorder=0)
            ax[i].grid(alpha=0.4)
            ax[i].plot(times[i], targ_flux_corr[i], '.', zorder=1, color=cmap(0), alpha=0.5)
            ax[i].errorbar(binned_times[i], binned_flux[i], binned_errs[i], linestyle='', marker='o', zorder=2, color=cmap(0))
            if i == 0:
                ax[i].set_ylabel('Flux', fontsize=14)
            else:
                ax[i].set_yticklabels([])
            ax[i].set_title(short_name, fontsize=20, color=cmap(0))
            ax[i].tick_params(labelsize=12)
            ax[i].set_ylim(1-standard_y, 1+standard_y)
            ax[i].set_xlim(np.mean(times[i])-standard_time_range/2, np.mean(times[i])+standard_time_range/2)
            ax[i].set_xlabel('Time (JD$_{UTC})$', fontsize=20)
        

        #Now, plot each reference star's corrected flux on every night. 
        color = cmap(95)

        for k in range(num_refs):
            ax = axis[k+1,:]
            for i in range(num_nights):
                weight = norm_night_weights[i][k]
                ax[i].axhline(1, color='k', alpha=0.7, lw=1, zorder=0)
                ax[i].grid(alpha=0.4)
                ax[i].plot(times[i], corr_flux[i][:,k], '.', zorder=1, color=color, alpha=0.5)
                ax[i].errorbar(binned_times[i], binned_ref_fluxes[i][:,k], binned_ref_flux_errs[i][:,k], linestyle='', marker='o', zorder=2, color=color)
                if i == 0:
                    ax[i].set_ylabel('Flux', fontsize=14)
                else:
                    ax[i].set_yticklabels([])
                ax[i].set_title('CS '+str(k+1)+', weight = {:1.2f}'.format(weight), fontsize=20, color=color)
                ax[i].tick_params(labelsize=12)
                ax[i].set_ylim(1-standard_y, 1+standard_y)
                ax[i].set_xlim(np.mean(times[i])-standard_time_range/2, np.mean(times[i])+standard_time_range/2)
                ax[i].set_xlabel('Time (JD$_{UTC})$', fontsize=20)

    elif num_nights == 1:
        #First, plot the target's corrected flux for each night. 
        ax = axis[0]
        ax.axhline(1, color='k', alpha=0.7, lw=1, zorder=0)
        ax.grid(alpha=0.4)
        ax.plot(times[0], targ_flux_corr[0], '.', zorder=1, color=cmap(0), alpha=0.5)
        ax.errorbar(binned_times[0], binned_flux[0], binned_errs[0], linestyle='', marker='o', zorder=2, color=cmap(0))
        ax.set_ylabel('Flux', fontsize=14)
        ax.set_title(short_name, fontsize=20, color=cmap(0))
        ax.tick_params(labelsize=12)
        ax.set_ylim(0.9,1.1)
        ax.set_xlim(np.mean(times[0])-standard_time_range/2, np.mean(times[0])+standard_time_range/2)
        ax.set_xlabel('Time (JD$_{UTC})$', fontsize=20)

        #Now, plot each reference star's corrected flux on every night. 
        color = cmap(95)
        for k in range(num_refs):
            ax = axis[k+1]
            weight = norm_night_weights[0][k]
            ax.axhline(1, color='k', alpha=0.7, lw=1, zorder=0)
            ax.grid(alpha=0.4)
            ax.plot(times[0], corr_flux[0][:,k], '.', zorder=1, color=color, alpha=0.5)
            ax.errorbar(binned_times[0], binned_ref_fluxes[0][:,k], binned_ref_flux_errs[0][:,k], linestyle='', marker='o', zorder=2, color=color)
            ax.set_ylabel('Flux', fontsize=14)
            ax.set_title('CS '+str(k+1)+', weight = {:1.2f}'.format(weight), fontsize=20, color=color)
            ax.tick_params(labelsize=12)
            ax.set_ylim(0.9,1.1)
            ax.set_xlim(np.mean(times[0])-standard_time_range/2, np.mean(times[0])+standard_time_range/2)
            ax.set_xlabel('Time (JD$_{UTC})$', fontsize=20)

    plt.tight_layout()
    
    #Save the figure. 
    print('Saving lightcurves for all sources...')
    if phot_type == 'aper':
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_'+'r='+ap_rad+'_nightly_corrected_flux_all_targets.png'
        if not os.path.exists(analysis_path/('aper_phot_analysis/'+ap_rad)):
            os.mkdir(analysis_path/('aper_phot_analysis/'+ap_rad))
        output_path = analysis_path/('aper_phot_analysis/'+ap_rad+'/'+filename)
        plt.savefig(output_path, dpi=300)
    else:
        raise RuntimeError("HAVE TO UPDATE!!!")
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_nightly_corrected_flux_all_targets.png'