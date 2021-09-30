from pines_analysis_toolkit.analysis.block_splitter import block_splitter
import numpy as np 

def block_binner(raw_times, raw_flux, time_threshold=0.15):
    block_inds = block_splitter(raw_times, time_threshold=time_threshold)
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