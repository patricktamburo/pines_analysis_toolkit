from astropy.io import fits
import numpy as np
import pdb
import glob 
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
import os
import time
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader
import getpass 
import paramiko
import pandas
import natsort
from datetime import datetime
import pysftp
from pathlib import Path
from progressbar import ProgressBar
from pines_analysis_toolkit.utils.quick_plot import quick_plot as qp

'''Authors:
		Paul Dalba, Boston University, February 2017
		Patrick Tamburo, Boston University, June-September 2020
	Purpose:
        Creates a master dome flat field image for a given date and band, and uploads to the PINES calibrations folder. 
        NOTE: Only admins are able to upload these files!
	Inputs:
        sftp (pysftp.Connection): the sftp connection to the pines server.
		date (str): the UT date during which the dome flat field data was obtained (i.e., '20200531')
        band (str): the band in which the flat field data was takend, 'J' or 'H'
        flat_start (int, optional): The file number that represents the start of the lights **on** dome flat sequence. Can use these arguments
            to specify the lights on/off data, in case the appropriate comments weren't entered into the FITS headers when they were created.
        flat_end (int, optional): The file number that represents the end of the lights **off** dome flat sequence.
            NOTE: This will only work if the lights_on sequence and lights_off sequence were taken back-to-back (which should always be the case).
        upload (bool, optional): Whether or not to upload the master dome flat to pines.bu.edu. By default, set to False (so you will not attempt to upload!).
        delete_raw (bool, optional): Whether or not to delete raw data from your local machine when the master flat is created. By default, set to False (won't delete raw dome flat images by default!).
    Outputs:
		None
	TODO:
		Tabular printing
        Bad image flagging
'''


def dome_flat_field(date, band, lights_on_start=0, lights_on_stop=0, lights_off_start=0, lights_off_stop=0, upload=False, delete_raw=False, sftp=''):
    clip_lvl = 3 #The value to use for sigma clipping. 
    np.seterr(invalid='ignore') #Suppress some warnings we don't care about in median combining. 
    plt.ion() #Turn on interactive plotting.
    pines_path = pines_dir_check()

    t1 = time.time()
    #If an sftp connection to the PINES server was passed, download the flat data. 
    if type(sftp) == pysftp.Connection:
        sftp.chdir('/')
        sftp.chdir('data/raw/mimir')
        run_list = sftp.listdir()
        data_path = '' #Initialize to check that it gets filled. 
        for i in range(len(run_list)):
            sftp.chdir(run_list[i])    
            date_list = sftp.listdir()
            if date in date_list:
                data_path = sftp.getcwd()
                print('{} directory found in pines.bu.edu:{}/'.format(date,data_path))
                print('')
                sftp.chdir(date)
                break
            sftp.chdir('..')
        
        if data_path == '':
            print('ERROR: {} not found in any run on pines.bu.edu:data/raw/mimir/.'.format(date))
            return
        
        else:
            #If the file start/stop numbers are specfied, grab those files.
            if (lights_on_stop != 0) or (lights_off_stop != 0):
                files_in_dir = sftp.listdir()
                on_flat_filenums = np.arange(lights_on_start, lights_on_stop+1, step=1)
                off_flat_filenums = np.arange(lights_off_start, lights_off_stop+1, step=1)
                flat_files = []
                lights_on_files = []
                lights_off_files = []

                #Add the lights-on flats to the file list. 
                for i in range(len(on_flat_filenums)):
                    file_num = on_flat_filenums[i]
                    #Generate the filename. 
                    if file_num < 10:
                        file_name = date+'.00'+str(file_num)+'.fits'
                    elif (file_num >= 10) and (file_num < 100):
                        file_name = date+'.0'+str(file_num)+'.fits'
                    else:
                        file_name = date+'.'+str(file_num)+'.fits'
                    #Check if the file name is in the directory, and if so, append it to the list of flat files. 
                    if file_name in files_in_dir:
                        flat_files.append(file_name)
                        lights_on_files.append(file_name)
                    else:
                        print('{} not found in directory, skipping.'.format(file_name))        

                #Do the same for the lights-off files. 
                for i in range(len(off_flat_filenums)):
                    file_num = off_flat_filenums[i]
                    #Generate the filename. 
                    if file_num < 10:
                        file_name = date+'.00'+str(file_num)+'.fits'
                    elif (file_num >= 10) and (file_num < 100):
                        file_name = date+'.0'+str(file_num)+'.fits'
                    else:
                        file_name = date+'.'+str(file_num)+'.fits'
                    #Check if the file name is in the directory, and if so, append it to the list of flat files. 
                    if file_name in files_in_dir:
                        flat_files.append(file_name)
                        lights_off_files.append(file_name)
                    else:
                        print('{} not found in directory, skipping.'.format(file_name)) 
            #Otherwise, find the files automatically using the night's log. 
            else:
                log_path = pines_path/'Logs'
                #Check if you already have the log for this date, if not, download it. 
                if not (log_path/(date+'_log.txt')).exists():
                    print('Downloading {}_log.txt to {}'.format(date,log_path))
                    sftp.get(date+'_log.txt',log_path/(date+'_log.txt'))

                log = pines_log_reader(log_path/(date+'_log.txt'))
                
                #Identify flat files. 
                flat_inds = np.where((log['Target'] == 'Flat') & (log['Filename'] != 'test.fits') & (log['Filt.'] == band))[0]
                flat_files = natsort.natsorted(list(set(log['Filename'][flat_inds]))) #Set guarantees we only grab the unique files that have been identified as flats, in case the log bugged out. 

            print('Found {} flat files.'.format(len(flat_files)))
            print('')

            #Downoad data to the appropriate Calibrations/Flats/Domeflats/band/Raw/ directory. 
            dome_flat_raw_path = pines_path/('Calibrations/Flats/Domeflats/'+band+'/Raw')
            for j in range(len(flat_files)):
                if not (dome_flat_raw_path/flat_files[j]).exists():
                    sftp.get(flat_files[j],(dome_flat_raw_path/flat_files[j]))
                    print('Downloading {} to {}, {} of {}.'.format(flat_files[j],dome_flat_raw_path, j+1, len(flat_files)))
                else:
                    print('{} already in {}, skipping download.'.format(flat_files[j],dome_flat_raw_path))
            print('')

            if (lights_on_stop == 0) and (lights_off_stop == 0):
                #Find the lights-on and lights-off flat files. 
                lights_on_files = []
                lights_off_files = []
                for j in range(len(flat_files)):
                    header = fits.open(dome_flat_raw_path/flat_files[j])[0].header
                    if header['FILTNME2'] != band:
                        print('ERROR: {} taken in filter other than {}. Double check your date, try specifying start/stop file numbers, etc.'.format(flat_files[j], band))
                        return
                    if header['OBJECT'] == 'dome_lamp_on':
                        lights_on_files.append(flat_files[j])
                    elif header['OBJECT'] == 'dome_lamp_off':
                        lights_off_files.append(flat_files[j])
                    else:
                        print("ERROR: header['OBJECT'] for {} is not 'dome_lamp_on' or 'dome_lamp_off. Double check your date, try specifying start/stop file numbers, etc.".format(flat_files[j]))
    else:
        dome_flat_raw_path = pines_path/('Calibrations/Flats/Domeflats/'+band+'/Raw')
        flat_files = natsort.natsorted(list(Path(dome_flat_raw_path).rglob(date+'*.fits')))
        #Find the lights-on and lights-off flat files. 
        lights_on_files = []
        lights_off_files = []
        for j in range(len(flat_files)):
            header = fits.open(dome_flat_raw_path/flat_files[j])[0].header
            if header['FILTNME2'] != band:
                print('ERROR: {} taken in filter other than {}. Double check your date, try specifying start/stop file numbers, etc.'.format(flat_files[j], band))
                return
            if header['OBJECT'] == 'dome_lamp_on':
                lights_on_files.append(flat_files[j])
            elif header['OBJECT'] == 'dome_lamp_off':
                lights_off_files.append(flat_files[j])
            else:
                print("ERROR: header['OBJECT'] for {} is not 'dome_lamp_on' or 'dome_lamp_off. Double check your date, try specifying start/stop file numbers, etc.".format(flat_files[j]))

    if len(lights_on_files) == 0 or len(lights_off_files) == 0:
        raise RuntimeError('No raw lights on/off flat files found with date '+date+' in '+band+' band!')

    print('Found {} lights-on flat files.'.format(len(lights_on_files)))
    print('Found {} lights-off flat files.'.format(len(lights_off_files)))
    print('')
    time.sleep(1)

    #Make cube of the lights-on images.
    num_images = len(lights_on_files)
    print('Reading in ', num_images,' lights-on flat images.')
    flat_lights_on_cube_raw = np.zeros([len(lights_on_files),1024,1024]) #Declare datatype to match raw mimir data. 
    print('')
    print('Flat frame information')
    print('-------------------------------------------------')
    print('ID   Mean               Stddev         Max    Min')
    print('-------------------------------------------------')
    for j in range(len(lights_on_files)):
        image_data = fits.open(dome_flat_raw_path/lights_on_files[j])[0].data[0:1024,:]
        header = fits.open(dome_flat_raw_path/lights_on_files[j])[0].header
        if header['FILTNME2'] != band:
            print('ERROR: {} taken in filter other than {}. Double check your date, try specifying start/stop file numbers, etc.'.format(lights_on_files[j], band))
            return
        flat_lights_on_cube_raw[j,:,:] = image_data #This line trims off the top two rows of the image, which are overscan.
        print(str(j+1)+'    '+str(np.mean(image_data))+'    '+str(np.std(image_data))+'    '+ str(np.std(image_data))+'    '+str(np.amin(image_data)))

    #For each pixel, calculate the mean, median, and standard deviation "through the stack" of lights on flat images.
    lights_on_cube_shape = np.shape(flat_lights_on_cube_raw)
    master_flat_lights_on = np.zeros((lights_on_cube_shape[1], lights_on_cube_shape[2]), dtype='float32')
    master_flat_lights_on_stddev = np.zeros((lights_on_cube_shape[1], lights_on_cube_shape[2]), dtype='float32')
    print('')
    print('Combining the lights-on flats.')
    print('......')
    pbar = ProgressBar()
    for x in pbar(range(lights_on_cube_shape[1])):
        for y in range(lights_on_cube_shape[2]):
            through_stack = flat_lights_on_cube_raw[:,y,x]
            through_stack_median = np.nanmedian(through_stack)
            through_stack_stddev = np.nanstd(through_stack)

            #Flag values that are > clip_lvl-sigma discrepant from the median.
            good_inds = np.where((abs(through_stack - through_stack_median) / through_stack_stddev <= clip_lvl))[0]

            #Calculate the sigma-clipped mean and sigma-clipped stddev using good_inds. 
            s_c_mean = np.nanmean(through_stack[good_inds])
            s_c_stddev = np.nanstd(through_stack[good_inds])

            #Store the sigma-clipped mean as the master dark value for this pixel. 
            master_flat_lights_on[y,x] = s_c_mean
            master_flat_lights_on_stddev[y,x] = s_c_stddev    

    #Make cube of the lights-off images
    num_images = len(lights_off_files)
    print('Reading in ', num_images,' lights-off flat images.')
    flat_lights_off_cube_raw = np.zeros([len(lights_off_files),1024,1024]) 
    print('')
    print('Flat frame information')
    print('-------------------------------------------------')
    print('ID   Mean               Stddev         Max    Min')
    print('-------------------------------------------------')
    for j in range(len(lights_off_files)):
        image_data = fits.open(dome_flat_raw_path/lights_off_files[j])[0].data[0:1024,:]
        header = fits.open(dome_flat_raw_path/lights_off_files[j])[0].header
        if header['FILTNME2'] != band:
            print('ERROR: {} taken in filter other than {}. Double check your date, try specifying start/stop file numbers, etc.'.format(lights_off_files[j], band))
            return
        flat_lights_off_cube_raw[j,:,:] = image_data #This line trims off the top two rows of the image, which are overscan.
        print(str(j+1)+'    '+str(np.mean(image_data))+'    '+str(np.std(image_data))+'    '+ str(np.std(image_data))+'    '+str(np.amin(image_data)))
    time.sleep(1)

    #For each pixel, calculate the mean, median, and standard deviation "through the stack" of lights off flat images.
    lights_off_cube_shape = np.shape(flat_lights_off_cube_raw)
    master_flat_lights_off = np.zeros((lights_off_cube_shape[1], lights_off_cube_shape[2]), dtype='float32')
    master_flat_lights_off_stddev = np.zeros((lights_off_cube_shape[1], lights_off_cube_shape[2]), dtype='float32')
    print('')
    print('Combining the lights-off flats.')
    print('......')
    pbar = ProgressBar()
    for x in pbar(range(lights_off_cube_shape[1])):
        for y in range(lights_off_cube_shape[2]):
            through_stack = flat_lights_off_cube_raw[:,y,x]
            through_stack_median = np.nanmedian(through_stack)
            through_stack_stddev = np.nanstd(through_stack)

            #Flag values that are > clip_lvl-sigma discrepant from the median.
            good_inds = np.where((abs(through_stack - through_stack_median) / through_stack_stddev <= clip_lvl))[0]

            #Calculate the sigma-clipped mean and sigma-clipped stddev using good_inds. 
            s_c_mean = np.nanmean(through_stack[good_inds])
            s_c_stddev = np.nanstd(through_stack[good_inds])

            #Store the sigma-clipped mean as the master dark value for this pixel. 
            master_flat_lights_off[y,x] = s_c_mean
            master_flat_lights_off_stddev[y,x] = s_c_stddev 
    
    #Create the master flat
    master_flat = master_flat_lights_on - master_flat_lights_off
    master_flat_norm = np.nanmedian(master_flat)
    master_flat = master_flat / master_flat_norm

    #Create the master error flat
    master_flat_error = np.sqrt(master_flat_lights_on_stddev**2 + master_flat_lights_off_stddev**2)
    master_flat_error = master_flat_error / master_flat_norm
    
    #Ourput flat files. 
    output_filename = pines_path/('Calibrations/Flats/Domeflats/'+band+'/Master Flats/master_flat_'+band+'_'+date+'.fits')
    #Add some header keywords detailing the master_dark creation process. 
    hdu = fits.PrimaryHDU(master_flat)
    hdu.header['HIERARCH DATE CREATED'] = datetime.utcnow().strftime('%Y-%m-%d')+'T'+datetime.utcnow().strftime('%H:%M:%S')

    #Now save to a file on your local machine. 
    print('')
    print('Writing the file to {}'.format(output_filename))
    #Check to see if other files of this name exist
    if os.path.exists(output_filename):
        print('')
        print('WARNING: This will overwrite {}!'.format(output_filename))
        flat_check = input('Do you want to continue? y/n: ')
        if flat_check == 'y':
            hdu.writeto(output_filename,overwrite=True)
            print('Wrote to {}!'.format(output_filename))
        else:
            print('Not overwriting!')
    else:
        hdu.writeto(output_filename,overwrite=True)
    print('')

    
    if upload:
        print('Beginning upload process to pines.bu.edu...')
        print('Note, only PINES admins will be able to upload.')
        time.sleep(2)
        print('')
        sftp.chdir('/')
        sftp.chdir('data/calibrations/Flats/Domeflats/'+band)
        upload_name = 'master_flat_'+band+'_'+date+'.fits'
        if upload_name in sftp.listdir():
            print('WARNING: This will overwrite {} in pines.bu.edu:data/calibrations/Flats/Domeflats/{}/'.format(upload_name,band))
            upload_check = input('Do you want to continue? y/n: ')
            if upload_check == 'y':
                sftp.put(output_filename,upload_name)
                print('Uploaded to pines.bu.edu:data/calibrations/Flats/Domeflats/{}/ !'.format(band))
            else:
                print('Skipping upload!')
        else:
            sftp.put(output_filename,upload_name)
            print('Uploaded {} to pines.bu.edu:data/calibrations/Flats/Domeflats/{}/!'.format(upload_name, band))


    output_filename = pines_path/('Calibrations/Flats/Domeflats/'+band+'/Master Flats Stddev/master_flat_stddev_'+band+'_'+date+'.fits')
    #Add some header keywords detailing the master_dark creation process. 
    hdu = fits.PrimaryHDU(master_flat_error)
    hdu.header['HIERARCH DATE CREATED'] = datetime.utcnow().strftime('%Y-%m-%d')+'T'+datetime.utcnow().strftime('%H:%M:%S')

    #Now save to a file on your local machine. 
    print('')
    print('Writing the file to {}'.format(output_filename))
    #Check to see if other files of this name exist
    if os.path.exists(output_filename):
        print('')
        print('WARNING: This will overwrite {}!'.format(output_filename))
        flat_check = input('Do you want to continue? y/n: ')
        if flat_check == 'y':
            hdu.writeto(output_filename,overwrite=True)
            print('Wrote to {}!'.format(output_filename))
        else:
            print('Not overwriting!')
    else:
        hdu.writeto(output_filename,overwrite=True)
    print('')

    
    if upload:
        print('Beginning upload process to pines.bu.edu...')
        print('Note, only PINES admins will be able to upload.')
        time.sleep(2)
        print('')
        sftp.chdir('/')
        sftp.chdir('data/calibrations/Flats Stddev/Domeflats/'+band)
        upload_name = 'master_flat_'+band+'_'+date+'.fits'
        if upload_name in sftp.listdir():
            print('WARNING: This will overwrite {} in pines.bu.edu:data/calibrations/Flats/Domeflats/{}/'.format(upload_name,band))
            upload_check = input('Do you want to continue? y/n: ')
            if upload_check == 'y':
                sftp.put(output_filename,upload_name)
                print('Uploaded to pines.bu.edu:data/calibrations/Flats Stddev/Domeflats/{}/ !'.format(band))
            else:
                print('Skipping upload!')
        else:
            sftp.put(output_filename,upload_name)
            print('Uploaded {} to pines.bu.edu:data/calibrations/Flats Stddev/Domeflats/{}/!'.format(upload_name, band))
    print('')
    if delete_raw:
        files_to_delete = glob.glob(os.path.join(dome_flat_raw_path/'*.fits'))
        for j in range(len(files_to_delete)):
            os.remove(files_to_delete[j])

    print('dome_flat_field runtime: ', np.round((time.time()-t1)/60,1), ' minutes.')
    print('Done!')