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

'''Authors:
		Paul Dalba, Boston University, February 2017
		Patrick Tamburo, Boston University, June 2020
	Purpose:
        Creates a master dome flat field image for a given run and band, and uploads to the PINES calibrations folder. 
        NOTE: Only admins are able to upload these files!
	Inputs:
		date (str): the UT date during which the dome flat field data was obtained (i.e., '20200531')
        band (str): the band in which the flat field data was takend, 'J' or 'H'
        lights_on_start (int, optional): The file number that represents the start of the lights on dome flat sequence. Can use these arguments
            to specify the lights on/off data, in case the appropriate comments weren't entered into the FITS headers when they were created.
        lights_on_end (int, optional): The file number that represents the end of the lights on dome flat sequence.
        lights_off_start (int, optional): The file number that represents the start of the lights off dome flat sequence.
        lights_off_end (int, optional): The file number that represents the end of the lights off dome flat sequence.
        upload (bool, optional): Whether or not to upload the master dome flat to pines.bu.edu. By default, set to True (so you will attempt to upload!).
        delete_raw (bool, optional): Whether or not to delete raw data from your local machine when the master flat is created. By default, set to False (won't delete raw dome flat images by default!).
    Outputs:
		None
	TODO:
		Tabular printing.
        Bad image flagging.
'''


def dome_flat_field(date, band, lights_on_start=0, lights_on_stop=0, lights_off_start=0, lights_off_stop=0, upload=True, delete_raw=False):
    clip_lvl = 3 #The value to use for sigma clipping. 
    np.seterr(invalid='ignore') #Suppress some warnings we don't care about in median combining. 
    plt.ion() #Turn on interactive plotting.

    t1 = time.time()

    pines_path = pines_dir_check()


    #Prompt login: 
    username = input('Enter username: ')
    password = getpass.getpass('Enter password: ')

     #Open ssh connection and set up local/remote paths.
    ssh = paramiko.SSHClient()
    ssh.load_system_host_keys()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect('pines.bu.edu',username=username, password=password)
    username = ''
    password = ''
    sftp = ssh.open_sftp()
    sftp.chdir('data/raw/mimir')
    run_list = sftp.listdir()
    data_path = '' #Initialize to check that it gets filled. 
    for i in range(len(run_list)):
        sftp.chdir(run_list[i])    
        date_list = sftp.listdir()
        if date in date_list:
            data_path = sftp.getcwd()
            print('{} directory found in pines.bu.edu:{}/'.format(date,data_path))
            sftp.chdir(date)
            break
        sftp.chdir('..')
    
    if data_path == '':
        print('ERROR: specified date not found in any run on pines.bu.edu:data/raw/mimir/')
        return
    
    else:
        #Download the log and read it in. 
        sftp.get(date+'_log.txt',pines_path+'Logs/'+date+'_log.txt')
        log = pines_log_reader(pines_path+'Logs/'+date+'_log.txt')

        #Identify flat files. 
        flat_inds = np.where((log['Target'] == 'Flat') & (log['Filename'] != 'test.fits'))[0]
        flat_files = natsort.natsorted(list(set(log['Filename'][flat_inds]))) #Set guarantees we only grab the unique files that have been identified as flats, in case the log bugged out. 
        print('Found {} flat files.'.format(len(flat_files)))
        print('')

        #Downoad this data to the appropriate Calibrations/Flats/Domeflats/band/Raw/ directory. 
        for j in range(len(flat_files)):
            if not os.path.exists(pines_path+'Calibrations/Flats/Domeflats/'+band+'/Raw/'+flat_files[j]):
                sftp.get(flat_files[j],pines_path+'Calibrations/Flats/Domeflats/'+band+'/Raw/'+flat_files[j])
                print('Downloading {} to {}, {} of {}.'.format(flat_files[j],pines_path+'Calibrations/Flats/Domeflats/'+band+'/Raw/', j+1, len(flat_files)))
            else:
                print('{} already in {}, skipping download.'.format(flat_files[j],pines_path+'Calibrations/Flats/Domeflats/'+band+'/Raw/'))
        print('')

        #Find the lights-on and lights-off flat files. 
        lights_on_files = []
        lights_off_files = []
        for j in range(len(flat_files)):
            header = fits.open(pines_path+'Calibrations/Flats/Domeflats/'+band+'/Raw/'+flat_files[j])[0].header
            if header['FILTNME2'] != band:
                print('ERROR: {} taken in filter other than {}. Double check your date, try specifying start/stop file numbers, etc.'.format(flat_files[j], band))
                return

            if header['OBJECT'] == 'dome_lamp_on':
                lights_on_files.append(flat_files[j])
            elif header['OBJECT'] == 'dome_lamp_off':
                lights_off_files.append(flat_files[j])
            else:
                print("ERROR: header['OBJECT'] for {} is not 'dome_lamp_on' or 'dome_lamp_off. Double check your date, try specifying start/stop file numbers, etc.".format(flat_files[j]))

        print('Found {} lights-on flat files.'.format(len(lights_on_files)))
        print('Found {} lights-off flat files.'.format(len(lights_off_files)))
        print('')
        time.sleep(1)

        #Make cube of the lights-on images.
        num_images = len(lights_on_files)
        print('Reading in ', num_images,' lights-on flat images.')
        flat_lights_on_cube_raw = np.zeros([len(lights_on_files),1024,1024]) 
        print('')
        print('Flat frame information')
        print('-------------------------------------------------')
        print('ID   Mean               Stddev         Max    Min')
        print('-------------------------------------------------')
        for j in range(len(lights_on_files)):
            image_data = fits.open(pines_path+'Calibrations/Flats/Domeflats/'+band+'/Raw/'+lights_on_files[j])[0].data[0:1024,:]
            flat_lights_on_cube_raw[j,:,:] = image_data #This line trims off the top two rows of the image, which are overscan.
            print(str(j+1)+'    '+str(np.mean(image_data))+'    '+str(np.std(image_data))+'    '+ str(np.std(image_data))+'    '+str(np.amin(image_data)))
        time.sleep(1)

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
            image_data = fits.open(pines_path+'Calibrations/Flats/Domeflats/'+band+'/Raw/'+lights_off_files[j])[0].data[0:1024,:]
            flat_lights_off_cube_raw[j,:,:] = image_data #This line trims off the top two rows of the image, which are overscan.
            print(str(j+1)+'    '+str(np.mean(image_data))+'    '+str(np.std(image_data))+'    '+ str(np.std(image_data))+'    '+str(np.amin(image_data)))
        time.sleep(1)


        #Combine flats into one master flat, which is bias-corrected, dark rate corrected, and normalized.

        #First, combine the lights-on cube
        lights_on_cube_shape = np.shape(flat_lights_on_cube_raw)
        #Now sigma clip the master flat during combining if desired. The clip level is set at the top of the code, and the sigma clipping will run iteratively until the clip no
        #longer removes any pixels or until the number of images drops by 50%.
        
        #The sigma clipping comes from Paul Dalba's original photometry code, and I'm trusting that he did it correctly. 
        print('')
        print('Combining the lights-on flats')
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

        #Subtract the lights off image from the lights on image.
        flat_master = lights_on_flat_master - lights_off_flat_master

        #Normalize the master flat.
        flat_master = flat_master / np.nanmedian(flat_master)

        #Plots
        plt.figure()
        avg,med,std = sigma_clipped_stats(flat_master)
        plt.imshow(flat_master,origin='lower',vmin=med-3*std,vmax=med+3*std)
        plt.title('Master Flat')
        plt.show()

        np.seterr(invalid='warn') #Turn invalid warnings back on, in case it would permanently turn it off otherwise.


        output_filename = pines_path+'Calibrations/Flats/Domeflats/'+band+'/Master Flats/master_flat_'+band+'_'+date+'.fits'
        
        #Now save to a file on your local machine. 
        print('')
        print('Writing the file to master_flat.fits')
        #Check to see if other files of this name exist
        if os.path.exists(output_filename):
            print('')
            print('WARNING: This will overwrite {}!'.format(output_filename))
            flat_check = input('Do you want to continue? y/n: ')
            if flat_check == 'y':
                fits.PrimaryHDU(flat_master).writeto(output_filename,overwrite=True)
                print('Wrote to {}!'.format(output_filename))
            else:
                print('Not overwriting!')
        else:
            fits.PrimaryHDU(flat_master).writeto(output_filename,overwrite=True)
        print('')

        
        if upload:
            print('Beginning upload process to pines.bu.edu...')
            time.sleep(2)
            print('')
            sftp.chdir('..')
            sftp.chdir('..')
            sftp.chdir('..')
            sftp.chdir('..')
            sftp.chdir('calibrations/Flats/Domeflats/'+band)
            upload_name = 'master_flat_'+band+'_'+date+'.fits'
            if upload_name in sftp.listdir():
                print('WARNING: This will overwrite {} in pines.bu.edu:data/calibrations/Flats/Domeflats/{}/'.format(upload_name,band))
                upload_check = input('Do you want to continue? y/n: ')
                if upload_check == 'y':
                    sftp.put(output_filename,'master_flat_'+band+'_'+date+'.fits')
                    print('Uploaded to pines.bu.edu:data/calibrations/Flats/Domeflats/{}/ !'.format(band))
                else:
                    print('Skipping upload!')
            else:
                sftp.put(output_filename,'master_flat_'+band+'_'+date+'.fits')
                print('Uploaded {} to pines.bu.edu:data/calibrations/Flats/Domeflats/{}/'.format(upload_name, band))

        print('')
        if delete_raw:
            files_to_delete = glob.glob(pines_path+'Calibrations/Flats/Domeflats/'+band+'/Raw/*.fits')
            for j in range(len(files_to_delete)):
                os.remove(files_to_delete[j])

        print('dome_flat_field runtime: ', np.round((time.time()-t1)/60,1), ' minutes.')
        print('Done!')