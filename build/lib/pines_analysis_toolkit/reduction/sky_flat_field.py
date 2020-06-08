from astropy.io import fits
import numpy as np
import pdb
import glob 
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
import os
import time

path = 'P:\\Photometry\\PINES\\202001\\skyflats\\'
dark = fits.open('P:\\Photometry\\PINES\\202001\\darks\\0.1\\master_dark_0.1.fits')[0].data
band = 'H'
sky_flat_files = np.sort(glob.glob(path+band+'\\*.fits'))
clip_lvl = 3

#Check if output directory is there, if not, make it.
output_path = path+band+'\\master_flat\\'
if not os.path.exists(output_path):
    os.mkdir(output_path)

#Make cube of the sky flats
num_images = len(sky_flat_files)
print('Combining ', num_images,' flat images into one master sky flat image.')
time.sleep(2)
sky_flat_cube_raw = np.zeros([len(sky_flat_files),1024,1024]) 
print('')
print('Flat frame information')
print('-------------------------------------------------')
print('ID   Mean               Stddev         Max    Min')
print('-------------------------------------------------')
for i in range(len(sky_flat_files)):
    image_data = fits.open(sky_flat_files[i])[0].data[0:1024,:] - dark #Read in, subtracting off dark. 
    sky_flat_cube_raw[i,:,:] = image_data #This line trims off the top two rows of the image, which are overscan.
    print(str(i)+'    '+str(np.mean(image_data))+'    '+str(np.std(image_data))+'    '+ str(np.amax(image_data))+'    '+str(np.amin(image_data)))

#Allow user to manually remove any bad flats.
print('')
print('Do any of the lights-on flat files need to be removed?')
check = input('Enter the ID #s as a list [] or return to continue: ')
if check != '':
    check = eval(check)
    good_set = list(set(range(len(sky_flat_files)))-set(check))
    sky_flat_cube_raw = sky_flat_cube_raw[good_set,:,:]

#Do we have to subtract a dark???

sky_flat_cube_shape =  np.shape(sky_flat_cube_raw)
"""Combine flats into one master flat, which is bias-corrected, dark rate corrected, and normalized."""

#Now sigma clip the master flat during combining if desired. The clip level is set by the
#user up in the input and the sigma clipping will run iteratively until the clip no
#longer removes any pixels or until the number of images drops by 50%.
print('')
print('Combining the sky flats')
print('Sigma Clipping starting')
print('(Ignore "invalid value in greater" warning)')
print('......')
#Start while loop that will iterate until clipping is finished
nan_count, clip_count = 0., 0.

while clip_count < 3.:
        avg, stddev = np.nanmean(sky_flat_cube_raw,axis=0), np.nanstd(sky_flat_cube_raw,axis=0)
        #Stddev must not be zero - if so, all elements are same value and we can plug in
        # an infinitesimal number instead of zero
        stddev[np.where(stddev==0.)] = 1e-4
        #Input the nans in the elements that should be clipped
        for i in range(sky_flat_cube_shape[0]):
            sky_flat_cube_raw[i,:,:][abs(sky_flat_cube_raw[i,:,:] - avg)/stddev > clip_lvl] = np.nan
        #Up the counter and keep track of the number of nans
        clip_count +=1 
        if np.sum(np.isnan(sky_flat_cube_raw)*1.) == nan_count: 
            break
        nan_count = np.sum(np.isnan(sky_flat_cube_raw)*1.)

#Combine the clipped cube
sky_flat_master, std_sky_flat_master = np.nanmean(sky_flat_cube_raw,axis=0), np.nanstd(sky_flat_cube_raw,axis=0)

flat_master = sky_flat_master / np.nanmedian(sky_flat_master)

plt.figure()
avg,med,std = sigma_clipped_stats(flat_master)
plt.imshow(flat_master,origin='lower',vmin=med-3*std,vmax=med+3*std)
plt.title('Master Sky Flat')
plt.show()

# domeflat = fits.open('P:\\Photometry\\PINES\\202001\\domeflats\\H\\master_flat\\master_flat_H.fits')[0].data

# diff_flat = flat_master-domeflat
# plt.figure()
# avg,med,std = sigma_clipped_stats(diff_flat)
# plt.imshow(diff_flat,origin='lower',vmin=med-3*std,vmax=med+3*std)
# plt.title('Master Sky Flat - Master Dome Flat')
# plt.show()

#Try reducing an image and compare it to domeflat_reduced image
path = 'P:\\Photometry\\PINES\\202001\\2MASS 0916+1057\\'
domeflat = fits.open(path+'domeflats\\master_flat_H.fits')[0].data
dark = fits.open(path+'darks\\master_dark_15.fits')[0].data
skyflat = flat_master
raw = fits.open(path+'raw\\20200202.048.fits')[0].data[0:1024,:]

skyflat_red = (raw-dark)/skyflat
domeflat_red = (raw-dark)/domeflat

avg,med,std = sigma_clipped_stats(skyflat_red)
print('std skyflat reduced: ',std)
avg2,med2,std2 = sigma_clipped_stats(domeflat_red)
print('std domeflat reduced: ',std2)

plt.ion()
fig,ax = plt.subplots(nrows=1,ncols=2,figsize=(10,8),sharex=True,sharey=True)
ax[0].imshow(skyflat_red,origin='lower',vmin=med,vmax=med+3*std)
ax[0].set_title('Skyflat reduced')
ax[1].imshow(domeflat_red,origin='lower',vmin=med2,vmax=med2+3*std2)
ax[1].set_title('Domeflat reduced')
pdb.set_trace()
#Now save to a file

'''#Now write the master_flat to a file'''
print('')
print('Writing the file to master_flat.fits')
#Check to see if other files of this name exist
if np.size(glob.glob(output_path+'master_flat_'+band+'.fits')) > 0.:
    print('')
    print('This will overwrite current master_flat_'+band+'.fits file in ',output_path,'!')
    flat_check = input('Hit return to continue, or ctrl+c to escape.')
fits.PrimaryHDU(flat_master).writeto(output_path+'master_skyflat_'+band+'.fits',overwrite=True)
print('File saved!')
pdb.set_trace()