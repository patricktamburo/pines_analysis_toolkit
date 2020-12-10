import matplotlib.pyplot as plt 
import numpy as np 
import pdb
import julian 
import matplotlib.dates as mdates
import batman 

def pines_lc_plot(times, targ_flux_corr, binned_times, binned_flux, binned_errs, analysis_path, time_mode='JD'):
    '''Authors:
        Patrick Tamburo, Boston University, November 2020
	Purpose:
        Creates a standard plot of corrected target flux. 
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


    #If UT, convert Julian dates to datetime objects. 
    if time_mode == 'UT':
        times = np.array([julian.from_jd(times[i]) for i in range(len(times))])
        binned_times =  np.array([julian.from_jd(binned_times[i]) for i in range(len(binned_times))])

    #Plot the corrected target flux. 
    fig, ax = plt.subplots(1,1,figsize=(14,5))
    ax.axhline(1, color='r', lw=2, zorder=0)
    #ax.plot(t_plot, f_plot, color='r', lw=2, zorder=0, label='1.5 R$_\oplus$ transit model')
    ax.grid(alpha=0.4)
    ax.plot(times, targ_flux_corr, linestyle='', marker='o', zorder=1, color='darkgrey', ms=5, label='Unbinned data')
    ax.errorbar(binned_times, binned_flux, binned_errs, linestyle='', marker='o', zorder=2, color='k', capsize=2, ms=7, label='Binned data')
    five_sig = 5*np.mean(binned_errs)
    ax.axhline(1-five_sig, color='k', lw=2, linestyle='--', zorder=0, label='5-$\sigma$ threshold')
    ax.tick_params(labelsize=14)
    ax.set_ylabel('Normalized Flux', fontsize=16)

    handles, labels = plt.gca().get_legend_handles_labels()
    order=[0,2,1]
    ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], fontsize=14, loc='best', framealpha=0.7)
    ax.set_ylim(0.9,1.1)
    if time_mode == 'UT':
        myFmt = mdates.DateFormatter('%H:%M')
        ax.xaxis.set_major_formatter(myFmt)
        fig.autofmt_xdate()     
        ax.set_xlabel('Time (UT)', fontsize=16)
    else:
        ax.set_xlabel('Time (JD$_{UTC}$)', fontsize=16)

    ax.set_title('{}, UT {} {}, {}'.format(short_name, times[0].strftime('%b'), times[0].day, times[0].year), fontsize=16)
    ax.set_xlim(times[0], times[-1])
    plt.tight_layout()

    plt.savefig(analysis_path/'lc.png', dpi=300)
    pdb.set_trace()