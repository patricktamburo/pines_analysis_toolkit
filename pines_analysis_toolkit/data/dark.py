import glob
import numpy as np
from astropy.io import fits
import pdb
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader
import paramiko
import time
import getpass
import natsort
import os
from datetime import datetime

'''Authors:
		Paul Dalba, Boston University, February 2017
		Patrick Tamburo, Boston University, June 2020
	Purpose:
        Creates a master dark image for a given date and exposure time, and uploads to the PINES calibrations folder. 
        NOTE: Only admins are able to upload these files!
	Inputs:
		date (str): the UT date during which the dome flat field data was obtained (i.e., '20200531')
        exptime (float): the exposure time of the dark images in question, in seconds
        dark_start (int, optional): The file number that represents the start of the dark sequence for this exptime. Can use these arguments
            to specify the dark data, in case the appropriate comments weren't entered into the FITS headers when they were created.
        dark_end (int, optional): The file number that represents the end of the dark sequence for this exptime.
        upload (bool, optional): Whether or not to upload the master dark to pines.bu.edu. By default, set to False (so you will not attempt to upload!).
        delete_raw (bool, optional): Whether or not to delete raw data from your local machine when the master dark is created. By default, set to False (won't delete raw dark images by default!).
    Outputs:
		None
	TODO:
		Tabular printing
        Bad image flagging
        Implement specified start/stop arguments

'''

def dark(date, exptime, dark_start=0, dark_stop=0, upload=False, delete_raw=False):
    clip_lvl = 3 #The value to use for sigma clipping. 
    pines_path = pines_dir_check()
    np.seterr(invalid='ignore') #Suppress some warnings we don't care about in median combining. 
    exptime = float(exptime)
    plt.ion() #Turn on interactive plotting.
    
    #Prompt login: 
    print('')
    username = input('Enter username: ')
    password = getpass.getpass('Enter password: ')

    t1 = time.time()

    #Open ssh connection and set up local/remote paths.
    ssh = paramiko.SSHClient()
    ssh.load_system_host_keys()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect('pines.bu.edu',username=username, password=password)
    sftp = ssh.open_sftp()
    password = ''
    print('')

    sftp.chdir('data/raw/mimir')
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
        log_path = pines_path/'Logs'
        #Check if you already have the log for this date, if not, download it. 
        if not (log_path/(date+'_log.txt')).exists():
            print('Donwloading {}_log.txt to {}\n'.format(date,log_path))
            sftp.get(date+'_log.txt',log_path/(date+'_log.txt'))
        
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

        num_images = len(dark_files)

        print('Reading in ', num_images,' dark images.')
        dark_cube_raw = np.zeros([len(dark_files),1024,1024]) 
        print('')
        print('Dark frame information')
        print('-------------------------------------------------')
        print('ID   Mean               Stddev         Max    Min')
        print('-------------------------------------------------')
        for j in range(len(dark_files)):
            image_data = fits.open(dark_path/dark_files[j])[0].data[0:1024,:] #This line trims off the top two rows of the image, which are overscan.
            dark_cube_raw[j,:,:] = image_data 
            print(str(j+1)+'    '+str(np.mean(image_data))+'    '+str(np.std(image_data))+'    '+ str(np.amax(image_data))+'    '+str(np.amin(image_data)))


        cube_shape = np.shape(dark_cube_raw)

        print('')
        print('Combining the darks')
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

        #Force to type float32 to preserve the size of raw Mimir images
        dark_master = dark_master.astype('float32')
               
        plt.ion()
        avg,med,std = sigma_clipped_stats(dark_master)
        plt.imshow(dark_master,origin='lower',vmin=med-3*std,vmax=med+3*std)
        plt.title('Master Dark, exptime = {} s'.format(str(exptime)))
        plt.show()

        np.seterr(invalid='warn') #Turn invalid warnings back on, in case it would permanently turn it off otherwise.


        output_filename = pines_path/('Calibrations/Darks/Master Darks/master_dark_'+str(exptime)+'_s_'+date+'.fits')
        
        #Add some header keywords detailing the master_dark creation process. 
        hdu = fits.PrimaryHDU(dark_master)
        hdu.header['HIERARCH MASTER_DARK CREATOR'] = username
        hdu.header['HIERARCH DATE CREATED'] = datetime.utcnow().strftime('%Y-%m-%d')+'T'+datetime.utcnow().strftime('%H:%M:%S')
        username = ''

        #Now save to a file on your local machine. 
        print('')
        print('Writing the file to master_dark_'+str(exptime)+'_s_'+date+'.fits')

        #Check to see if other files of this name exist
        if os.path.exists(output_filename):
            print('')
            print('WARNING: This will overwrite {}!'.format(output_filename))
            dark_check = input('Do you want to continue? y/n: ')
            if dark_check == 'y':
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
            sftp.chdir('..')
            sftp.chdir('..')
            sftp.chdir('..')
            sftp.chdir('..')
            sftp.chdir('calibrations/Darks')
            upload_name = 'master_dark_'+str(exptime)+'_s_'+date+'.fits'
            if upload_name in sftp.listdir():
                print('WARNING: This will overwrite {} in pines.bu.edu:data/calibrations/Darks/'.format(upload_name))
                upload_check = input('Do you want to continue? y/n: ')
                if upload_check == 'y':
                    sftp.put(output_filename,upload_name)
                    print('Uploaded to pines.bu.edu:data/calibrations/Darks/!')
                else:
                    print('Skipping upload!')
            else:
                sftp.put(output_filename,upload_name)
                print('Uploaded {} to pines.bu.edu:data/calibrations/Darks/!'.format(upload_name))

        print('')
        if delete_raw:
            files_to_delete = glob.glob(os.path.join(dark_path/'*.fits'))
            for j in range(len(files_to_delete)):
                os.remove(files_to_delete[j])
        sftp.close()
        print('dark runtime: ', np.round((time.time()-t1)/60,1), ' minutes.')
        print('Done!')