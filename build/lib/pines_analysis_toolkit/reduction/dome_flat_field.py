from astropy.io import fits
import numpy as np
import pdb
import glob 
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
import os
import time

'''Purpose: 
        Makes master domeflat using lights-on and lights-off images.
   Assumed directory structure: 
        A top-level directory that contains directories for the bands you're using (i.e. J/ or H/, or both).
        The J/ and H/ directories should contain directories called lights-on/ and lights-off/, which contain the raw lights-on and lights-off flat images for that band.
        The program will save a file called master_flat_band.fits, where band is J or H. This file gets saved to a directory called master_flat/, which gets created at
        the same level as the lights-on/ and lights-off directories. 
'''


#----------------------------------------------USER INPUTS----------------------------------------------------
band = 'H' #The band you're working in. 
run = '202003' #The run the flat data is from. 
#-------------------------------------------------------------------------------------------------------------

path = '/Users/tamburo/Documents/Data/PINES/Calibrations/Flats/Domeflats/'+band


lights_on_flat_files = np.sort(glob.glob(path+'/Raw data/'+run+'/lights-on/*.fits'))
lights_off_flat_files = np.sort(glob.glob(path+'/Raw data/'+run+'/lights-off/*.fits'))
clip_lvl = 3

#Check if output directory is there, if not, make it.
output_path = path+'/Master flats/' #Saves 
if not os.path.exists(output_path):
    os.mkdir(output_path)

#Make cube of the lights-on images
num_images = len(lights_on_flat_files)
print('Combining ', num_images,' flat images into one master lights-on flat image.')
time.sleep(2)
flat_lights_on_cube_raw = np.zeros([len(lights_on_flat_files),1024,1024]) 
print('')
print('Flat frame information')
print('-------------------------------------------------')
print('ID   Mean               Stddev         Max    Min')
print('-------------------------------------------------')
for i in range(len(lights_on_flat_files)):
    image_data = fits.open(lights_on_flat_files[i])[0].data[0:1024,:]
    flat_lights_on_cube_raw[i,:,:] = image_data #This line trims off the top two rows of the image, which are overscan.
    print(str(i)+'    '+str(np.mean(image_data))+'    '+str(np.std(image_data))+'    '+ str(np.amax(image_data))+'    '+str(np.amin(image_data)))

#Allow user to manually remove any bad flats.
print('')
print('Do any of the lights-on flat files need to be removed?')
check = input('Enter the ID #s as a list [] or return to continue: ')
if check != '':
    check = eval(check)
    good_set = list(set(range(len(lights_on_flat_files)))-set(check))
    flat_lights_on_cube_raw = flat_lights_on_cube_raw[good_set,:,:]

#Make cube of the lights-off images
num_images = len(lights_off_flat_files)
print('Combining ', num_images,' flat images into one master lights-off flat image.')
time.sleep(2)
flat_lights_off_cube_raw = np.zeros([len(lights_off_flat_files),1024,1024]) 
print('')
print('Flat frame information')
print('-------------------------------------------------')
print('ID   Mean               Stddev         Max    Min')
print('-------------------------------------------------')
for i in range(len(lights_off_flat_files)):
    image_data = fits.open(lights_off_flat_files[i])[0].data[0:1024,:]
    flat_lights_off_cube_raw[i,:,:] = image_data #This line trims off the top two rows of the image, which are overscan.
    print(str(i)+'    '+str(np.mean(image_data))+'    '+str(np.std(image_data))+'    '+ str(np.amax(image_data))+'    '+str(np.amin(image_data)))

#Allow user to manually remove any bad flats.
print('')
print('Do any of the lights-on flat files need to be removed?')
check = input('Enter the ID #s as a list [] or return to continue: ')
if check != '':
    check = eval(check)
    good_set = list(set(range(len(lights_off_flat_files)))-set(check))
    flat_lights_off_cube_raw = flat_lights_off_cube_raw[good_set,:,:]
    

    pdb.set_trace()
"""Combine flats into one master flat, which is bias-corrected, dark rate corrected, and normalized."""

#First, combine the lights-on cube
lights_on_cube_shape = np.shape(flat_lights_on_cube_raw)
#Now sigma clip the master flat during combining if desired. The clip level is set by the
#user up in the input and the sigma clipping will run iteratively until the clip no
#longer removes any pixels or until the number of images drops by 50%.
print('')
print('Combining the lights-on flats')
print('Sigma Clipping starting')
print('(Ignore "invalid value in greater" warning)')
print('......')
#Start while loop that will iterate until clipping is finished
nan_count, clip_count = 0., 0.

while clip_count < 3.:
        avg, stddev = np.nanmean(flat_lights_on_cube_raw,axis=0), np.nanstd(flat_lights_on_cube_raw,axis=0)
        #Stddev must not be zero - if so, all elements are same value and we can plug in
        # an infinitesimal number instead of zero
        stddev[np.where(stddev==0.)] = 1e-4
        #Input the nans in the elements that should be clipped
        for i in range(lights_on_cube_shape[0]):
            flat_lights_on_cube_raw[i,:,:][abs(flat_lights_on_cube_raw[i,:,:] - avg)/stddev > clip_lvl] = np.nan
        #Up the counter and keep track of the number of nans
        clip_count +=1 
        if np.sum(np.isnan(flat_lights_on_cube_raw)*1.) == nan_count: 
            break
        nan_count = np.sum(np.isnan(flat_lights_on_cube_raw)*1.)

#Combine the clipped cube
lights_on_flat_master, std_lights_on_flat_master = np.nanmean(flat_lights_on_cube_raw,axis=0), np.nanstd(flat_lights_on_cube_raw,axis=0)

#Then, combine the lights-off cube
lights_off_cube_shape = np.shape(flat_lights_off_cube_raw)
#Now sigma clip the master flat during combining if desired. The clip level is set by the
#user up in the input and the sigma clipping will run iteratively until the clip no
#longer removes any pixels or until the number of images drops by 50%.
print('')
print('Combining the lights-off flats')
print('Sigma Clipping starting')
print('(Ignore "invalid value in greater" warning)')
print('......')
#Start while loop that will iterate until clipping is finished
nan_count, clip_count = 0., 0.

while clip_count < 3.:
        avg, stddev = np.nanmean(flat_lights_off_cube_raw,axis=0), np.nanstd(flat_lights_off_cube_raw,axis=0)
        #Stddev must not be zero - if so, all elements are same value and we can plug in
        # an infinitesimal number instead of zero
        stddev[np.where(stddev==0.)] = 1e-4
        #Input the nans in the elements that should be clipped
        for i in range(lights_off_cube_shape[0]):
            flat_lights_off_cube_raw[i,:,:][abs(flat_lights_off_cube_raw[i,:,:] - avg)/stddev > clip_lvl] = np.nan
        #Up the counter and keep track of the number of nans
        clip_count +=1 
        if np.sum(np.isnan(flat_lights_off_cube_raw)*1.) == nan_count: 
            break
        nan_count = np.sum(np.isnan(flat_lights_off_cube_raw)*1.)
#Combine the clipped cube
lights_off_flat_master, std_lights_off_flat_master = np.nanmean(flat_lights_off_cube_raw,axis=0), np.nanstd(flat_lights_off_cube_raw,axis=0)

# plt.ion()
# plt.figure()
# avg,med,std = sigma_clipped_stats(lights_on_flat_master)
# plt.imshow(lights_on_flat_master,origin='lower',vmin=med-3*std,vmax=med+3*std)
# plt.title('Lights-on Image')
# plt.show()

# plt.figure()
# avg,med,std = sigma_clipped_stats(lights_off_flat_master)
# plt.imshow(lights_off_flat_master,origin='lower',vmin=med-3*std,vmax=med+3*std)
# plt.title('Lights-off Image')
# plt.show()

#lights_on_flat_master = np.nanmedian(flat_lights_on_cube_raw,axis=0)
#lights_off_flat_master = np.nanmedian(flat_lights_off_cube_raw,axis=0)
flat_master = lights_on_flat_master - lights_off_flat_master
flat_master = flat_master / np.nanmedian(flat_master)

# plt.figure()
# avg,med,std = sigma_clipped_stats(flat_master)
# plt.imshow(flat_master,origin='lower',vmin=med-3*std,vmax=med+3*std)
# plt.title('Master Flat')
# plt.show()

#Get the ut date when the first image was taken, and append it to the filename. 
ut_date_str = fits.open(lights_on_flat_files[0])[0].header['DATE-OBS'].split('T')[0].replace('-','') 

#Now save to a file
'''#Now write the master_flat to a file'''
print('')
print('Writing the file to master_flat.fits')
#Check to see if other files of this name exist
if np.size(glob.glob(output_path+'master_flat_'+band+'_'+ut_date_str+'.fits')) > 0.:
    print('')
    print('This will overwrite current master_flat_'+band+'.fits file in ',output_path,'!')
    flat_check = input('Hit return to continue, or ctrl+c to escape.')
fits.PrimaryHDU(flat_master).writeto(output_path+'master_flat_'+band+'_'+ut_date_str+'.fits',overwrite=True)
print('File saved!')
