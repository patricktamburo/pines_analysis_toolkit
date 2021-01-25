import numpy as np
import pdb

'''Authors:
		Patrick Tamburo, Boston University, June 2020
	Purpose:
        Finds individual blocks of data given a single night of exposture times. 
	Inputs:
        times (numpy array): 1D array of exposure times (in days)
        bad_vals (numpy array, optional): 1D array of any sigma-clipped values to be excluded in the averaging/error estimation 
    Outputs:
        night_inds (list): list containing the data indices associated with each night
	TODO:
'''
def block_splitter(times, bad_vals=[]):
    block_boundaries = np.where(np.gradient(times) > 0.15/24)[0]

    #Staring observations. TODO: This will not work if there is a mix of staring/hopping observations on a single night!
    if len(block_boundaries) == 0:
        time_bin = 10.0 #Minutes over which to bin. 
        block_boundaries = np.where(np.gradient(((times - times[0]) * (24*60)) % time_bin) < 0.2)[0]
    
    num_blocks = int(1 + len(block_boundaries) / 2)
    block_inds = [[] for x in range(num_blocks)]
    #Hopping observations.
    for j in range(num_blocks):
        if j == 0:
            block_inds[j].extend(np.arange(0,block_boundaries[0]+1))
        elif j > 0 and j < num_blocks - 1:
            block_inds[j].extend(np.arange(block_boundaries[2*j-1],block_boundaries[2*j]+1))
        else:
            block_inds[j].extend(np.arange(block_boundaries[2*j-1],len(times)))
    if len(bad_vals) > 0:
        for i in range(len(block_inds)):
            bad_locs = []
            for j in range(len(block_inds[i])):
                if block_inds[i][j] in bad_vals:
                    bad_locs.append(np.where(block_inds[i] == block_inds[i][j])[0][0])
            if len(bad_locs) != 0:
                bad_locs = bad_locs[::-1]
                for j in range(len(bad_locs)):
                    block_inds[i].remove(block_inds[i][bad_locs[j]])
    return block_inds
