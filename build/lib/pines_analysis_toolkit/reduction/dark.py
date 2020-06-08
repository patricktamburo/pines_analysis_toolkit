import glob
import numpy as np
from astropy.io import fits
import pdb
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats

path = 'P:\\Photometry\\PINES\\202001\\darks\\60\\'
exptime = '60'
dark_files = np.sort(glob.glob(path+'*.fits'))
num_images = len(dark_files)
clip_lvl = 3

print('Combining ', num_images,' dark images into one master dark image.')
dark_cube_raw = np.zeros([len(dark_files),1024,1024]) 
print('')
print('Dark frame information')
print('-------------------------------------------------')
print('ID   Mean               Stddev         Max    Min')
print('-------------------------------------------------')
for i in range(len(dark_files)):
    image_data = fits.open(dark_files[i])[0].data[0:1024,:] #This line trims off the top two rows of the image, which are overscan.
    dark_cube_raw[i,:,:] = image_data 
    print(str(i)+'    '+str(np.mean(image_data))+'    '+str(np.std(image_data))+'    '+ str(np.amax(image_data))+'    '+str(np.amin(image_data)))

#Allow user to manually remove any bad darkes
print('')
print('Do any of the dark files need to be removed?')
check = input('Enter the ID #s as a list [] or return to continue: ')
if check != '':
    check = eval(check)
    good_set = list(set(range(len(dark_files)))-set(check))
    dark_cube_raw = dark_cube_raw[good_set,:,:]

cube_shape = np.shape(dark_cube_raw)

#Now sigma clip the master dark during combining if desired. The clip level is set by the
# user up in the input and the sigma clipping will run iteratively until the clip no
# longer removes any pixels or until the number of images drops by 50%.
print('')
print('Sigma Clipping starting')
print('(Ignore "invalid value in greater" warning)')
print('......')
#Start while loop that will iterate until clipping is finished
nan_count, clip_count = 0., 0.

while clip_count < 3.:
    avg, stddev = np.nanmean(dark_cube_raw,axis=0), np.nanstd(dark_cube_raw,axis=0)
    #Stddev must not be zero - if so, all elements are same value and we can plug in
    # an infinitesimal number instead of zero
    stddev[np.where(stddev==0.)] = 1e-4
    #Input the nans in the elements that should be clipped
    for i in range(cube_shape[0]):
        dark_cube_raw[i,:,:][abs(dark_cube_raw[i,:,:] - avg)/stddev > clip_lvl] = np.nan
    #Up the counter and keep track of the number of nans
    clip_count +=1 
    if np.sum(np.isnan(dark_cube_raw)*1.) == nan_count: 
        break
    nan_count = np.sum(np.isnan(dark_cube_raw)*1.)

#Combine the clipped cube
dark_master, std_dark_master = np.nanmean(dark_cube_raw,axis=0), np.nanstd(dark_cube_raw,axis=0)
plt.ion()
avg,med,std = sigma_clipped_stats(dark_master)
plt.imshow(dark_master,origin='lower',vmin=med-3*std,vmax=med+3*std)
plt.show()
pdb.set_trace()
'''#Now write the master_flat to a file'''
print('')
print('Writing the file to master_dark_'+exptime+'.fits')
#Check to see if other files of this name exist
if np.size(glob.glob(path+'master_dark_'+exptime+'.fits')) > 0.:
    print('')
    print('This will overwrite current master_dark_'+exptime+'.fits file in ',path,'!')
    dark_check = input('Hit return to continue, or ctrl+c to escape.')
fits.PrimaryHDU(dark_master).writeto(path+'master_dark_'+exptime+'.fits',overwrite=True)
print('File saved!')
pdb.set_trace()
