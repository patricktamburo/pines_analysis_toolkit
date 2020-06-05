import pdb
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import glob
import pickle
import os

def bias():
    data_path = 'P:\\Data\\201911\\'
    date_list = os.listdir(data_path)
    clip_level = 3

    #Get list of bias files in these directories. 
    bias_list = []
    for date in date_list:
        path = data_path+date+'\\'
        files = np.sort(glob.glob(path+'*.fits'))
        for file in files:
            header = fits.open(file)[0].header
            if (header['OBJECT'] == 'BIAS') and (header['EXPTIME'] == 0):
                bias_list.append(file)
        print('Done ',date)
    print('Found ',len(bias_list),' bias files.')

    #Declare an empty cube into which bias data will be read. 
    bias_cube_raw = np.zeros([len(bias_list),1024,1024]) 

    for i in range(len(bias_list)):
        image_data = fits.open(bias_list[i])[0].data
        bias_cube_raw[i,:,:] = image_data[0:1024,:]
        print(str(i)+'    '+str(np.mean(image_data))+'    '+str(np.std(image_data))+'    '+ str(np.amax(image_data))+'    '+str(np.amin(image_data)))

    #Allow user to manually remove any bad biases
    print('')
    print('Do any of the bias files need to be removed?')
    check = input('Enter the ID #s as a list [] or return to continue: ')
    if check != '':
        check = eval(check)
        good_set = list(set(range(len(bias_list)))-set(check))
        bias_cube_raw = bias_cube_raw[good_set,:,:]

    #Now sigma clip the master bias during combining if desired. The clip level is set by the
    # user up in the input and the sigma clipping will run iteratively until the clip no
    # longer removes any pixels or until the number of images drops by 50%.
    print('')
    print('Sigma Clipping starting')
    print('(Ignore "invalid value in greater" warning)')
    print('......')
    #Start while loop that will iterate until clipping is finished
    nan_count, clip_count = 0., 0.

    while clip_count < 3.:
        avg, stddev = np.nanmean(bias_cube_raw,axis=0), np.nanstd(bias_cube_raw,axis=0)
        #Stddev must not be zero - if so, all elements are same value and we can plug in
        # an infinitesimal number instead of zero
        stddev[np.where(stddev==0.)] = 1e-4
        #Input the nans in the elements that should be clipped
        for i in range(len(bias_cube_raw)):
            bias_cube_raw[i,:,:][abs(bias_cube_raw[i,:,:] - avg)/stddev > clip_level] = np.nan
        #Up the counter and keep track of the number of nans
        clip_count +=1 
        if np.sum(np.isnan(bias_cube_raw)*1.) == nan_count: 
            break
        nan_count = np.sum(np.isnan(bias_cube_raw)*1.)
    #Combine the clipped cube
    bias_master, std_bias_master = np.nanmean(bias_cube_raw,axis=0), np.nanstd(bias_cube_raw,axis=0)
    print('Saving bias master to ',data_path,'master_bias.fits.')
    hdu = fits.PrimaryHDU(bias_master)
    hdu.writeto(data_path+'master_bias.fits')
    pdb.set_trace()
