from astropy.io import fits
import numpy as np
import pdb
import pandas as pd
import itertools as it

#Reads in two bias files, subtracts them, finds the standard deviation of the difference image, then measures read noise. 
#Following http://spiff.rit.edu/classes/phys445/lectures/readout/readout.html.
def readnoise_calculator():
    G = 8.21 #e- / DN

    bias_numbers = np.arange(23,43,1) #biases from this run were images 23 - 42 on 2019/11/14. 
    unique_combos = pd.Series(list(it.combinations(np.unique(bias_numbers),2))) #Get all unique combinations of those images. 
    read_noises = np.zeros(len(unique_combos))
    for i in range(len(unique_combos)):
        bias1_path = 'P:\\Data\\201911\\20191114\\20191114.0'+str(unique_combos[i][0])+'.fits'
        bias2_path = 'P:\\Data\\201911\\20191114\\20191114.0'+str(unique_combos[i][1])+'.fits'
        bias1 = fits.open(bias1_path)[0].data
        bias2 = fits.open(bias2_path)[0].data
        diff_image = bias1-bias2 #Make the difference image. 
        read_noises[i] = np.std(diff_image)/np.sqrt(2) * G #Calculate RN, converted to e-. 
        print(read_noises[i])
    pdb.set_trace()
