import numpy as np
import pdb 

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