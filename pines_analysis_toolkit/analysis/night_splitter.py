import numpy as np
import pdb 

'''Authors:
		Patrick Tamburo, Boston University, June 2020
	Purpose:
        Finds different nights of data given a 1D array of exposure times (in days), and returns their indices. 
	Inputs:
        times (numpy array): 1D array of exposure times (in days)
    Outputs:
        night_inds (list): list containing the data indices associated with each night
	TODO:
'''
def night_splitter(times):
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