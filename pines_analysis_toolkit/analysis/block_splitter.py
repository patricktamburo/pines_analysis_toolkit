import numpy as np
import pdb

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
