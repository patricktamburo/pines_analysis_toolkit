import matplotlib.pyplot as plt 
import numpy as np 
import pdb
import julian 
import matplotlib.dates as mdates

plt.ioff() 

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

    fig, ax = plt.subplots(nrows=1, ncols=num_nights, figsize=(17,5), sharey=True)
    for i in range(num_nights):
        ax[i].set_yscale('log')
        for j in range(num_refs):
            color = cm(1.*j/num_refs) 
            marker = markers[j % len(markers)]
            ax[i].errorbar(times[i], raw_ref_flux[i][:,j], raw_ref_err[i][:,j], marker=marker, linestyle='-', label='Ref. '+str(j+1), color=color)
        if i == 0:
            ax[i].set_ylabel('Photons', fontsize=16)
        ax[i].tick_params(labelsize=14)
        ax[i].set_xlabel('Time (JD$_{UTC}$)', fontsize=16)

        ax[i].errorbar(times[i], raw_targ_flux[i], raw_targ_err[i], marker='o', linestyle='-', label='Target', mfc='none', mew=2, color='tab:orange')

        if i == num_nights - 1:
            ax[i].legend(bbox_to_anchor=(1.1, 1), fontsize=14)
        
    plt.suptitle(short_name+' Nightly Raw Flux', fontsize=16)
    plt.subplots_adjust(left=0.07, wspace=0.1, top=0.93, bottom=0.17)

    if phot_type == 'aper':
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_'+'r='+ap_rad+'_nightly_raw_flux.png'
    else:
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_nightly_raw_flux.png'
    output_path = analysis_path/filename
    plt.savefig(output_path, dpi=300)

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
    ax.set_ylabel('Photons', fontsize=16)
    ax.tick_params(labelsize=14)
    ax.set_xlabel('Time (JD$_{UTC}$)', fontsize=16)

    ax.errorbar(times, raw_targ_flux[0], raw_targ_err[0], marker='o', linestyle='-', label='Target', mfc='none', mew=2, color='tab:orange')

    ax.legend(bbox_to_anchor=(1.1, 1), fontsize=14)
        
    plt.suptitle(short_name+' Global Raw Flux', fontsize=16)
    plt.subplots_adjust(left=0.07, wspace=0.1, top=0.93, bottom=0.17)

    if phot_type == 'aper':
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_'+'r='+ap_rad+'_global_raw_flux.png'
    else:
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_global_raw_flux.png'
    output_path = analysis_path/filename
    plt.savefig(output_path, dpi=300)

def normalized_flux_plot(times, norm_targ_flux, norm_targ_err, norm_ref_flux, norm_ref_err, alc_flux, short_name, analysis_path, phot_type, ap_rad):
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

    cm = plt.get_cmap('viridis_r')
    markers = ['+', 'x', '*', 'X']

    fig, ax = plt.subplots(nrows=1, ncols=num_nights, figsize=(17,5), sharey=True)
    for i in range(num_nights):
        for j in range(num_refs):
            color = cm(1.*j/num_refs) 
            marker = markers[j % len(markers)]
            ax[i].errorbar(times[i], norm_ref_flux[i][:,j], norm_ref_err[i][:,j], marker=marker, linestyle='', label='Ref. '+str(j+1), color=color)
        if i == 0:
            ax[i].set_ylabel('Normalized Flux', fontsize=16)
        ax[i].tick_params(labelsize=14)
        ax[i].set_xlabel('Time (JD$_{UTC}$)', fontsize=16)

        ax[i].errorbar(times[i], norm_targ_flux[i], norm_targ_err[i], marker='o', linestyle='-', label='Target', mfc='none', mew=2, color='tab:orange')
        ax[i].plot(times[i], alc_flux[i], color='tab:blue', lw=2, marker='o', label='ALC')
        if i == num_nights - 1:
            ax[i].legend(bbox_to_anchor=(1.05, 1), fontsize=14)
        
    plt.suptitle(short_name+' Nightly Normalized Flux', fontsize=16)
    plt.subplots_adjust(left=0.07, wspace=0.15, top=0.93, bottom=0.17)

    if phot_type == 'aper':
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_'+'r='+ap_rad+'_nightly_normalized_flux.png'
    else:
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_nightly_normalized_flux.png'
    output_path = analysis_path/filename
    plt.savefig(output_path, dpi=300)

def global_normalized_flux_plot(times, norm_targ_flux, norm_targ_err, norm_ref_flux, norm_ref_err, alc_flux, short_name, analysis_path, phot_type, ap_rad):
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

    cm = plt.get_cmap('viridis_r')
    markers = ['+', 'x', '*', 'X']

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(17,5))
    for j in range(num_refs):
        color = cm(1.*j/num_refs) 
        marker = markers[j % len(markers)]
        ax.errorbar(times, norm_ref_flux[0][:,j], norm_ref_err[0][:,j], marker=marker, linestyle='', label='Ref. '+str(j+1), color=color)
    ax.set_ylabel('Normalized Flux', fontsize=16)
    ax.tick_params(labelsize=14)
    ax.set_xlabel('Time (JD$_{UTC}$)', fontsize=16)

    ax.errorbar(times, norm_targ_flux[0], norm_targ_err[0], marker='o', linestyle='-', label='Target', mfc='none', mew=2, color='tab:orange')
    ax.plot(times, alc_flux[0], color='tab:blue', lw=2, marker='o', label='ALC')
    ax.legend(bbox_to_anchor=(1.1, 1), fontsize=14)
        
    plt.suptitle(short_name+' Global Normalized Flux', fontsize=16)
    plt.subplots_adjust(left=0.07, wspace=0.15, top=0.93, bottom=0.17)

    if phot_type == 'aper':
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_'+'r='+ap_rad+'_global_normalized_flux.png'
    else:
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_global_normalized_flux.png'
    output_path = analysis_path/filename
    plt.savefig(output_path, dpi=300)

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

    #Plot the corrected target flux. 
    fig, ax = plt.subplots(nrows=1,ncols=num_nights,figsize=(17,5), sharey=True)
    plt.subplots_adjust(left=0.07, wspace=0.1, top=0.86, bottom=0.17, right=0.96)

    for i in range(num_nights):
        ax[i].axhline(1, color='r', lw=2, zorder=0)
        ax[i].grid(alpha=0.4)
        ax[i].plot(times[i], targ_flux_corr[i], linestyle='', marker='o', zorder=1, color='darkgrey', ms=5, label='Unbinned data')
        ax[i].errorbar(binned_times[i], binned_flux[i], binned_errs[i], linestyle='', marker='o', zorder=2, color='k', capsize=2, ms=7, label='Binned data')
        five_sig = 5*np.mean(binned_errs[i])
        ax[i].axhline(1-five_sig, color='k', lw=2, linestyle='--', zorder=0, label='5-$\sigma$ threshold')
        ax[i].tick_params(labelsize=14)
        
        if i == 0:
            ax[i].set_ylabel('Normalized Flux', fontsize=14)

        # if i == num_nights - 1:
        #     ax[i].legend(bbox_to_anchor=(1.1, 1), fontsize=14)

        ax[i].set_xlabel('Time (JD$_{UTC}$)', fontsize=14)
        ut_date = julian.from_jd(times[i][0])
        ut_date_str = 'UT '+ut_date.strftime('%b. %d %Y')
        ax[i].set_title(ut_date_str, fontsize=14)

    plt.suptitle(short_name+' Nightly Corrected Target Lightcurve', fontsize=16)

    if phot_type == 'aper':
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_'+'r='+ap_rad+'_nightly_corrected_target_flux.png'
    else:
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_nightly_corrected_target_flux.png'
    plt.savefig(analysis_path/filename, dpi=300)


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
    plt.subplots_adjust(left=0.07, wspace=0.1, top=0.86, bottom=0.17, right=0.96)

    ax.axhline(1, color='r', lw=2, zorder=0)
    ax.grid(alpha=0.4)
    ax.plot(times, targ_flux_corr[0], linestyle='', marker='o', zorder=1, color='darkgrey', ms=5, label='Unbinned data')
    ax.errorbar(binned_times, binned_flux[0], binned_errs[0], linestyle='', marker='o', zorder=2, color='k', capsize=2, ms=7, label='Binned data')
    five_sig = 5*np.mean(binned_errs[0])
    ax.axhline(1-five_sig, color='k', lw=2, linestyle='--', zorder=0, label='5-$\sigma$ threshold')
    ax.tick_params(labelsize=14)
    
    ax.set_ylabel('Normalized Flux', fontsize=14)

    ax.set_xlabel('Time (JD$_{UTC}$)', fontsize=14)
    plt.suptitle(short_name+' Global Corrected Target Lightcurve', fontsize=16)

    if phot_type == 'aper':
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_'+'r='+ap_rad+'_global_corrected_target_flux.png'
    else:
        filename = short_name.replace(' ','')+'_'+phot_type+'_phot_global_corrected_target_flux.png'
    plt.savefig(analysis_path/filename, dpi=300)


def corr_all_sources_plot(times, targ_flux_corr, binned_times, binned_flux, binned_errs, corr_flux, binned_ref_fluxes, binned_ref_flux_errs, short_name, analysis_path, phot_type, ap_rad):
    pdb.set_trace()
    cmap = plt.cm.viridis
    fig, ax = plt.subplots(num_refs+1, 1, figsize=(14,5*num_refs))
    plt.subplots_adjust(top=0.95, left=0.05, bottom=0.05, right=0.95)

    ax[0].axhline(1, color='k', alpha=0.7, lw=1, zorder=0)
    ax[0].grid(alpha=0.4)
    ax[0].plot(times, targ_flux_corr, '.', zorder=1, color=cmap(0), alpha=0.5)
    ax[0].errorbar(binned_times, binned_flux, binned_errs, linestyle='', marker='o', zorder=2, color=cmap(0))
    ax[0].set_ylabel('Corrected Flux', fontsize=14)
    ax[0].set_title(short_name, fontsize=16, color=cmap(0))
    ax[0].tick_params(labelsize=12)
    ax[0].set_ylim(0.9,1.1)
    ax[0].set_xlabel('Time (JD$_{UTC})$', fontsize=14)
    
    for k in range(num_refs):
        color_ind = int(((k+2) * 256 / (num_refs+1))) - 1
        color = cmap(color_ind)
        ax[k+1].axhline(1, color='k', alpha=0.7, lw=1, zorder=0)
        ax[k+1].grid(alpha=0.4)
        ax[k+1].plot(times, corr_flux[:,k], '.', zorder=1, color=color, alpha=0.5)
        ax[k+1].errorbar(binned_times, binned_ref_fluxes[:,k], binned_ref_flux_errs[:,k], linestyle='', marker='o', zorder=2, color=color)
        ax[k+1].set_ylabel('Corrected Flux', fontsize=14)
        ax[k+1].set_title('CS '+str(k+1), fontsize=16, color=color)
        ax[k+1].tick_params(labelsize=12)
        ax[k+1].set_ylim(0.9,1.1)
        ax[k+1].set_xlabel('Time (JD$_{UTC})$', fontsize=14)
    plt.tight_layout()
    plt.savefig(analysis_path/'cs_fig.png', dpi=300)