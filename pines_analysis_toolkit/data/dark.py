import glob
import numpy as np
from astropy.io import fits
import pdb
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader
from pines_analysis_toolkit.utils import quick_plot as qps
import paramiko
import time
import getpass
import natsort
import os
from datetime import datetime
import pysftp
from pathlib import Path
from progressbar import ProgressBar

'''Authors:
		Paul Dalba, Boston University, February 2017
		Patrick Tamburo, Boston University, June-September 2020
	Purpose:
        Creates a master dark image for a given date and exposure time, and uploads to the PINES calibrations folder. 
        NOTE: Only admins are able to upload these files!
	Inputs:
		sftp (pysftp.Connection): the sftp connection to the pines server. 
        date (str): the UT date during which the dome flat field data was obtained (i.e., '20200531')
        exptime (float): the exposure time of the dark images in question, in seconds
        dark_start (int, optional): The file number that represents the start of the dark sequence for this exptime. Can use these arguments
            to specify the dark data, in case the appropriate comments weren't entered into the FITS headers when they were created.
        dark_end (int, optional): The file number that represents the end of the dark sequence for this exptime.
        upload (bool, optional): Whether or not to upload the master dark to pines.bu.edu. By default, set to False (so you will not attempt to upload!).
        delete_raw (bool, optional): Whether or not to delete raw data from your local machine when the master dark is created. By default, set to False (won't delete raw dark images by default!).
    Outputs:
		Saves master_dark_exptime_s_date.fits to Calibrations/Master Darks/
        Saves master_dark_stddev_exptime_s_date.fits to Calibrations/Master Darks Stddev/
	TODO:
        None
'''

def dark(date, exptime, dark_start=0, dark_stop=0, upload=False, delete_raw=False, sftp=''):
    clip_lvl = 3 #The value to use for sigma clipping. 
    pines_path = pines_dir_check()
    np.seterr(invalid='ignore') #Suppress some warnings we don't care about in median combining. 
    exptime = float(exptime)
    plt.ion() #Turn on interactive plotting.

    t1 = time.time()
    #If an sftp connection to the PINES server was passed, download the dark data. 
    if type(sftp) == pysftp.Connection:
        sftp.chdir('/data/raw/mimir')
        run_list = sftp.listdir()
        data_path = '' #Initialize to check that it gets filled. 
        for i in range(len(run_list)):
            sftp.chdir(run_list[i])    
            date_list = sftp.listdir()
            if date in date_list:
                data_path = sftp.getcwd()
                print('{} directory found in pines.bu.edu:{}/\n'.format(date,data_path))
                sftp.chdir(date)
                break
            sftp.chdir('..')
    
        if data_path == '':
            print('ERROR: specified date not found in any run on pines.bu.edu:data/raw/mimir/\n')
            return
        else:
            #If the file start/stop numbers are specified, grab those files.
            if (dark_stop != 0):
                files_in_dir = sftp.listdir()
                dark_filenums = np.arange(dark_start, dark_stop+1, step=1)
                dark_files = []

                #Add the darks to the file list. 
                for i in range(len(dark_filenums)):
                    file_num = dark_filenums[i]
                    #Generate the filename. 
                    if file_num < 10:
                        file_name = date+'.00'+str(file_num)+'.fits'
                    elif (file_num >= 10) and (file_num < 100):
                        file_name = date+'.0'+str(file_num)+'.fits'
                    else:
                        file_name = date+'.'+str(file_num)+'.fits'
                    #Check if the file name is in the directory, and if so, append it to the list of flat files. 
                    if file_name in files_in_dir:
                        dark_files.append(file_name)
                    else:
                        print('{} not found in directory, skipping.'.format(file_name))        
                pdb.set_trace()
            else:
                #Otherwise, find the files automatically using the night's log. 
                log_path = pines_path/'Logs'
                #Check if you already have the log for this date, if not, download it. 
                #Download from the /data/logs/ directory on PINES.
                if not (log_path/(date+'_log.txt')).exists():
                    print('Downloading {}_log.txt to {}\n'.format(date,log_path))
                    sftp.get('/data/logs/'+date+'_log.txt',log_path/(date+'_log.txt'))
                
                #Read in the log from this date.
                log = pines_log_reader(log_path/(date+'_log.txt'))

                #Identify dark files. 
                dark_inds = np.where((log['Target'] == 'Dark') & (log['Filename'] != 'test.fits') & (log['Exptime'] == exptime))[0]
                dark_files = natsort.natsorted(list(set(log['Filename'][dark_inds]))) #Set guarantees we only grab the unique files that have been identified as flats, in case the log bugged out. 
            print('Found {} dark files.'.format(len(dark_files)))
            print('')
            #Downoad data to the Calibrations/Darks/Raw/ directory. 
            dark_path = pines_path/('Calibrations/Darks/Raw')
            for j in range(len(dark_files)):
                if not (dark_path/dark_files[j]).exists():
                    sftp.get(dark_files[j],os.path.join(pines_path,dark_path/dark_files[j]))
                    print('Downloading {} to {}, {} of {}.'.format(dark_files[j], dark_path, j+1, len(dark_files)))
                else:
                    print('{} already in {}, skipping download.'.format(dark_files[j],dark_path))
            print('')

    #If no sftp was passed, search for files on disk. 
    else:
        dark_path = pines_path/('Calibrations/Darks/Raw')
        all_dark_files = natsort.natsorted(list(Path(dark_path).rglob(date+'*.fits')))
        dark_files = []
        for file in all_dark_files:
            if fits.open(file)[0].header['EXPTIME'] == exptime:
                dark_files.append(file)

    num_images = len(dark_files)

    if num_images == 0:
        raise RuntimeError('No raw dark files found on disk with date '+date+'!')

    print('Reading in ', num_images,' dark images.')
    dark_cube_raw = np.zeros([len(dark_files),1024,1024]) 
    print('')
    print('Dark frame information')
    print('-------------------------------------------------')
    print('ID   Mean               Stddev         Max    Min')
    print('-------------------------------------------------')
    for j in range(len(dark_files)):
        image_data = fits.open(dark_path/dark_files[j])[0].data[0:1024,:] #This line trims off the top two rows of the image, which are overscan.
        header = fits.open(dark_path/dark_files[j])[0].header
        if header['EXPTIME'] != exptime:
            print('ERROR: {} taken has exposure time different than than exptime.'.format(dark_files[j]))
            return
        dark_cube_raw[j,:,:] = image_data 
        print(str(j+1)+'    '+str(np.mean(image_data))+'    '+str(np.std(image_data))+'    '+ str(np.amax(image_data))+'    '+str(np.amin(image_data)))

    cube_shape = np.shape(dark_cube_raw)

    master_dark = np.zeros((cube_shape[1], cube_shape[2]), dtype='float32')
    master_dark_stddev = np.zeros((cube_shape[1], cube_shape[2]), dtype='float32')

    print('')
    print('Combining the darks')
    print('......')

    pbar = ProgressBar()
    #For each pixel, calculate the mean, median, and standard deviation "through the stack" of darks.
    for x in pbar(range(cube_shape[1])):
        for y in range(cube_shape[2]):
            through_stack = dark_cube_raw[:,y,x]
            through_stack_median = np.nanmedian(through_stack)
            through_stack_stddev = np.nanstd(through_stack)

            #Flag values that are > clip_lvl-sigma discrepant from the median.
            good_inds = np.where((abs(through_stack - through_stack_median) / through_stack_stddev <= clip_lvl))[0]

            #Calculate the sigma-clipped mean and sigma-clipped stddev using good_inds. 
            s_c_mean = np.nanmean(through_stack[good_inds])
            s_c_stddev = np.nanstd(through_stack[good_inds])

            #Store the sigma-clipped mean as the master dark value for this pixel. 
            master_dark[y,x] = s_c_mean
            master_dark_stddev[y,x] = s_c_stddev
    

    np.seterr(invalid='warn') #Turn invalid warnings back on, in case it would permanently turn it off otherwise.

    output_filename = pines_path/('Calibrations/Darks/Master Darks/master_dark_'+str(exptime)+'_s_'+date+'.fits')
    
    #Add some header keywords detailing the master_dark creation process. 
    hdu = fits.PrimaryHDU(master_dark)
    hdu.header['HIERARCH DATE CREATED'] = datetime.utcnow().strftime('%Y-%m-%d')+'T'+datetime.utcnow().strftime('%H:%M:%S')

    #Now save to a file on your local machine. 
    #Check to see if other files of this name exist.
    # if os.path.exists(output_filename):
    #     print('')
    #     print('WARNING: This will overwrite {}!'.format(output_filename))
    #     dark_check = input('Do you want to continue? y/n: ')
    #     if dark_check == 'y':
    #         hdu.writeto(output_filename,overwrite=True)
    #         print('Wrote to {}!'.format(output_filename))
    #     else:
    #         print('Not overwriting!')
    # else:
    hdu.writeto(output_filename,overwrite=True)
    print('Wrote to {}!'.format(output_filename))
    
    #Upload the master dark to PINES server.
    if upload:
        print('Uploading to pines.bu.edu...')
        sftp.chdir('..')
        sftp.chdir('..')
        sftp.chdir('..')
        sftp.chdir('..')
        sftp.chdir('calibrations/Darks')
        upload_name = 'master_dark_'+str(exptime)+'_s_'+date+'.fits'
        # if upload_name in sftp.listdir():
        #     print('WARNING: This will overwrite {} in pines.bu.edu:data/calibrations/Darks/'.format(upload_name))
        #     upload_check = input('Do you want to continue? y/n: ')
        #     if upload_check == 'y':
        #         sftp.put(output_filename,upload_name)
        #         print('Uploaded to pines.bu.edu:data/calibrations/Darks/!')
        #     else:
        #         print('Skipping upload!')
        # else:
        sftp.put(output_filename,upload_name)
        print('Uploaded {} to pines.bu.edu:data/calibrations/Darks/!'.format(upload_name))

        sftp.chdir('..')

    #Do the same thing for the sigma-clipped standard deviation image. 

    output_filename = pines_path/('Calibrations/Darks/Master Darks Stddev/master_dark_stddev_'+str(exptime)+'_s_'+date+'.fits')
    
    #Add some header keywords detailing the master_dark creation process. 
    hdu = fits.PrimaryHDU(master_dark_stddev)
    hdu.header['HIERARCH DATE CREATED'] = datetime.utcnow().strftime('%Y-%m-%d')+'T'+datetime.utcnow().strftime('%H:%M:%S')
    username = ''

    #Now save to a file on your local machine. 
    print('')

    #Check to see if other files of this name exist.
    # if os.path.exists(output_filename):
    #     print('')
    #     print('WARNING: This will overwrite {}!'.format(output_filename))
    #     dark_check = input('Do you want to continue? y/n: ')
    #     if dark_check == 'y':
    #         hdu.writeto(output_filename,overwrite=True)
    #         print('Wrote to {}!'.format(output_filename))
    #     else:
    #         print('Not overwriting!')
    # else:
    hdu.writeto(output_filename,overwrite=True)
    print('Wrote to {}!'.format(output_filename))
    
    #Upload the master dark to PINES server.
    if upload:
        print('Uploading to pines.bu.edu...')
        sftp.chdir('Darks Stddev')
        upload_name = 'master_dark_stddev_'+str(exptime)+'_s_'+date+'.fits'
        # if upload_name in sftp.listdir():
        #     print('WARNING: This will overwrite {} in pines.bu.edu:data/calibrations/Darks Stddev/'.format(upload_name))
        #     upload_check = input('Do you want to continue? y/n: ')
        #     if upload_check == 'y':
        #         sftp.put(output_filename,upload_name)
        #         print('Uploaded to pines.bu.edu:data/calibrations/Darks Stddev/!')
        #     else:
        #         print('Skipping upload!')
        # else:
        sftp.put(output_filename,upload_name)
        print('Uploaded {} to pines.bu.edu:data/calibrations/Darks Stddev/!'.format(upload_name))

    print('')
    #Delete raw dark images from disk.
    if delete_raw:
        files_to_delete = glob.glob(os.path.join(dark_path/'*.fits'))
        for j in range(len(files_to_delete)):
            os.remove(files_to_delete[j])

    print('dark runtime: ', np.round((time.time()-t1)/60,1), ' minutes.')
    print('Done!')