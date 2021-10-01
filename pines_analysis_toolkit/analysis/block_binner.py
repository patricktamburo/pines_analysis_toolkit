from pines_analysis_toolkit.analysis.block_splitter import block_splitter
import numpy as np 

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