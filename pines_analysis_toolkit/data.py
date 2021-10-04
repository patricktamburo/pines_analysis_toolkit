import pines_analysis_toolkit as pat 
from pines_analysis_toolkit.utils import pines_dir_check, pines_log_reader, short_name_creator, object_directory_creator, quick_plot as qp, pines_login

from astropy.stats import SigmaClip
from astropy.io import fits
from astropy.stats import histogram, sigma_clipped_stats

from scipy.stats import sigmaclip

from photutils import Background2D, MedianBackground
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import glob 
import time
import os
import pysftp
from pathlib import Path
from progressbar import ProgressBar
from natsort import natsorted

import matplotlib.pyplot as plt 

def bg_2d(data, box_size=50):
    """Removes large-scale changes in the background of an image. See https://photutils.readthedocs.io/en/stable/background.html.

    :param data: 2D array of pixel values
    :type data: numpy array
    :param box_size: Size of boxes used to do background estimation, defaults to 50
    :type box_size: int, optional
    :return: 2D array of pixel values with the background model subtracted
    :rtype: numpy array
    """
    sigma_clip = SigmaClip(sigma=3.)
    bkg_estimator = MedianBackground()
    bkg = Background2D(data, (box_size, box_size), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)

    return data - bkg.background

def bpm_chooser(bpm_path, header):
    """Identifies the closest bpm in time to the image you want reduce, that has the correct exposure time and is in the right band.

    :param bpm_path: path to the ~/PINES_analysis_toolkit/Calibrations/Bad Pixel Masks/ directory.
    :type bpm_path: pathlib.PosixPath
    :param header: header of the image you want to reduce
    :type header: astropy.io fits header
    :return: path to and name of chosen bad pixel mask 
    :rtype: pathlib.PosixPath, str
    """
    exptime = header['EXPTIME']
    band = header['FILTNME2']
    obs_date = datetime.strptime(header['DATE-OBS'].split('T')[0].replace('-',''),'%Y%m%d')
    possible_bpms = [x for x in (bpm_path.glob('*'+band+'_'+str(exptime)+'*.fits'))]

    if (len(possible_bpms) == 0): 
        print('Warning: Could not find any suitable bpms to reduce {}.'.format(header['FILENAME']))
        print('Expanding bpm search.')
        possible_bpms = [x for x in (bpm_path.glob('*'+band+'*.fits'))]
        if len(possible_bpms) == 0:
            print('ERROR: no suitable BPMs found inspect manually.')
            return
    possible_bpm_dates = [datetime.strptime(i.name.split('_')[-1].split('.')[0],'%Y%m%d') for i in possible_bpms]
    bpm_date_distances = [abs(possible_bpm_dates[i]-obs_date) for i in range(len(possible_bpm_dates))]
    bpm_ind = np.where(np.array(bpm_date_distances) == min(np.array(bpm_date_distances)))[0][0]
    master_bpm = fits.open(possible_bpms[bpm_ind])[0].data	
    master_bpm_name = possible_bpms[bpm_ind].name
    return master_bpm, master_bpm_name

def bpm_maker(flat_date, dark_date, exptime, band, upload=False, sftp=''):
    """Creates a combined bad pixel mask from Kokopelli, variable, hot, and dead pixel masks.

    :param flat_date: date of the master flat used to reduce the image
    :type flat_date: str
    :param dark_date: date of the master dark used to reduce the image
    :type dark_date: str
    :param exptime: exposure time in seconds
    :type exptime: float
    :param band: band of observations (e.g., 'J' or 'H')
    :type band: str
    :param upload: whether or not to upload the resulting bpm to the PINES server, defaults to False
    :type upload: bool, optional
    :param sftp: sftp connection to the PINES server, only needed if upload == True, defaults to ''
    :type sftp: str, optional
    """
    pines_path = pines_dir_check()

    #Load in the different masks. 
    kokopelli_path = pines_path/('Calibrations/Kokopelli Mask/kokopelli_mask.fits')
    kokopelli_mask = (1-fits.open(kokopelli_path)[0].data).astype('int')[0:1024,:]

    variable_path = pines_path/('Calibrations/Variable Pixel Masks/vpm_'+str(exptime)+'_s_'+dark_date+'.fits')
    variable_mask = fits.open(variable_path)[0].data

    hot_path = pines_path/('Calibrations/Hot Pixel Masks/hpm_'+str(exptime)+'_s_'+dark_date+'.fits')
    hot_mask = fits.open(hot_path)[0].data

    dead_path = pines_path/('Calibrations/Dead Pixel Masks/dpm_'+band+'_'+flat_date+'.fits')
    dead_mask = fits.open(dead_path)[0].data

    #Visualize all masks
    bpm = np.zeros(np.shape(dead_mask), dtype='int')
    bad_locs = np.where((kokopelli_mask == 1) | (variable_mask == 1) | (hot_mask == 1) | (dead_mask == 1))
    bpm[bad_locs] = 1
    
    num_bad = len(np.where(bpm ==1)[0])
    frac_bad = num_bad / 1024**2

    print('{} percent of the detector flagged as bad.'.format(np.round(frac_bad*100,1)))

    output_filename = 'bpm_'+band+'_'+str(exptime)+'_s_'+flat_date+'.fits'
    output_path = pines_path/('Calibrations/Bad Pixel Masks/'+output_filename)

    hdu = fits.PrimaryHDU(bpm)
    hdu.header['HIERARCH DATE CREATED'] = datetime.utcnow().strftime('%Y-%m-%d')+'T'+datetime.utcnow().strftime('%H:%M:%S')

    #Now save to a file on your local machine. 
    print('')
    print('Writing the file to '+output_filename)


    hdu.writeto(output_path,overwrite=True)
    print('Wrote to {}!'.format(output_path))
    print('')

    #Upload the master dark to PINES server.
    if upload:
        print('Beginning upload process to pines.bu.edu...')
        print('Note, only PINES admins will be able to upload.')
        print('')
        sftp.chdir('/')
        sftp.chdir('data/calibrations/Bad Pixel Masks')
        upload_name = output_filename
        sftp.put(output_path,upload_name)
        print('Uploaded {} to pines.bu.edu:data/calibrations/Bad Pixel Masks/!'.format(upload_name))
   
def dark(date, exptime, dark_start=0, dark_stop=0, upload=False, delete_raw=False, sftp=''):
    """ Creates a master dark image for a given date and exposure time, and uploads to the PINES calibrations folder

    :param date: date on which dark files were taken (YYYYMMDD)
    :type date: str
    :param exptime: exposure time in seconds
    :type exptime: str
    :param dark_start: file number indicating the start of a sequence of darks you want to use, defaults to 0
    :type dark_start: int, optional
    :param dark_stop: file number indicating the start of a sequence of darks you want to use, defaults to 0
    :type dark_stop: int, optional
    :param upload: whether or not to upload to the PINES server, defaults to False
    :type upload: bool, optional
    :param delete_raw: whether or not to delete raw dark files, defaults to False
    :type delete_raw: bool, optional
    :param sftp: sftp connection to the PINES server, defaults to ''
    :type sftp: str, optional
    :raises RuntimeError: breaks if no dark files are found
    """
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
                dark_files = natsorted(list(set(log['Filename'][dark_inds]))) #Set guarantees we only grab the unique files that have been identified as flats, in case the log bugged out. 
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
        all_dark_files = natsorted(list(Path(dark_path).rglob(date+'*.fits')))
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

def dead_pixels(date, band, upload=False, sftp=''):
    """Creates a dead pixel mask from a master flat image. 

    :param date: date of the master flat (YYYYMMDD)
    :type date: str
    :param band: band of flat image (e.g. 'J' or 'H')
    :type band: str
    :param upload: whether or not to upload the resulting dead pixel mask to the PINES server, defaults to False
    :type upload: bool, optional
    :param sftp: sftp connection to the PINES server, defaults to ''
    :type sftp: str, optional
    """
    pines_path = pines_dir_check()
    flat_path = pines_path/('Calibrations/Flats/Domeflats/'+band+'/Master Flats/')
    flat_file = natsorted(list(flat_path.rglob('*'+date+'.fits')))[0]
    master_flat = fits.open(flat_file)[0].data
    shape = np.shape(master_flat)

    # #Can plot a histogram of master_flat pixel values. 
    # his = histogram(master_flat, bins=250)
    # plt.figure(figsize=(12,5))
    # plt.bar(his[1][0:-1], his[0], width=0.005) 
    # plt.yscale('log')
    # plt.xlabel('Pixel Value', fontsize=14)
    # plt.ylabel('N$_{pix}$', fontsize=14)
    # plt.title('Distribution of Pixel Values in Master Flat '+band+' '+date, fontsize=16)
    # pdb.set_trace()

    #Incorporate bad pixel information from the Kokopelli pixel mask. 
    kokopelli_path = pines_path/('Calibrations/Kokopelli Mask/kokopelli_mask.fits')
    kokopelli_mask = (1-fits.open(kokopelli_path)[0].data).astype('int')[0:1024,:]
    master_flat[np.where(kokopelli_mask == 1)] = np.nan
    
    print('')
    print('Flagging dead pixels.')
    
    #Find hot pixels in the master dark. This uses the 8 pixels surrounding each pixel, ignoring the pixels on the edges of the detector.
    #This will iterate until no new dead pixels are found. 
    dead_pixel_mask = np.zeros((shape[0], shape[1]), dtype='int')
    total_flagged = 0

    box_ls = [13, 11, 9, 7, 5]
    clip_lvls = [4.5, 4.5, 4.5, 4, 4]
    iteration = 1
    for i in range(len(box_ls)):
        box_l = box_ls[i]
        clip_lvl = clip_lvls[i]
        print('Box size = {} x {}.'.format(box_l, box_l))
        print('Sigma clipping level = {}.'.format(clip_lvl))
        print('......')
        num_flagged = 999 #Initialize
        box_minus = int(box_l/2)
        box_plus = int(box_l/2) + 1
        targ_pix_ind = int(box_l**2/2)

        while (num_flagged > 1): 
            num_flagged = 0
            pbar = ProgressBar()
            for xx in pbar(range(int(box_l/2), shape[0]-int(box_l/2))):
                for yy in range(int(box_l/2), shape[1]-int(box_l/2)):
                    if (not np.isnan(master_flat[yy,xx])): #Only check pixels that aren't already flagged. 
                        box_2d = master_flat[yy-box_minus:yy+box_plus,xx-box_minus:xx+box_plus] #2d cutout surrounding the target pixel.
                        box_1d = box_2d[~np.isnan(box_2d)].ravel() #Unraveled box, ignoring any NaNs. 
                        neighbor_vals = np.delete(box_1d, targ_pix_ind) #Remove the target pixel from the 1d array. 

                        #neighbor_vals = sigmaclip(n, low=3.5, high=3.5)[0] #Clip away other bad pixels in the box to avoid biasing the mean/standard deviation.                   
                        if master_flat[yy,xx] < (np.mean(neighbor_vals) - clip_lvl*np.std(neighbor_vals)): #Flag pixel as hot if it is more than clip_lvl sigma higher than the mean of its neighbors.
                            dead_pixel_mask[yy,xx] = 1
                            master_flat[yy,xx] = np.nan #Set this pixel to a NaN so that it's ignored on subsequent iterations. 
                            num_flagged += 1
            iteration += 1
            total_flagged += num_flagged
            print('Iteration {}: {} new dead pixels identified, {} dead pixels total.'.format(iteration, num_flagged, total_flagged))
            print('')

    print('')
    print('Found {} dead pixels.'.format(total_flagged))

    output_filename = 'dpm_'+band+'_'+date+'.fits'
    output_path = pines_path/('Calibrations/Dead Pixel Masks/'+output_filename)

    hdu = fits.PrimaryHDU(dead_pixel_mask)
    hdu.header['HIERARCH DATE CREATED'] = datetime.utcnow().strftime('%Y-%m-%d')+'T'+datetime.utcnow().strftime('%H:%M:%S')
    hdu.header['HIERARCH SIGMA CLIP LVL'] = clip_lvl

    #Now save to a file on your local machine. 
    print('')
    print('Writing the file to '+output_filename)

    #Check to see if other files of this name exist.
    # if os.path.exists(output_path):
    #     print('')
    #     print('WARNING: This will overwrite {}!'.format(output_path))
    #     dark_check = input('Do you want to continue? y/n: ')
    #     if dark_check == 'y':
    #         hdu.writeto(output_path,overwrite=True)
    #         print('Wrote to {}!'.format(output_path))
    #     else:
    #         print('Not overwriting!')
    # else:
    hdu.writeto(output_path,overwrite=True)
    print('Wrote to {}!'.format(output_path))
    print('')

    #Upload the master dark to PINES server.
    if upload:
        print('Beginning upload process to pines.bu.edu...')
        print('Note, only PINES admins will be able to upload.')
        print('')
        sftp.chdir('/')
        sftp.chdir('data/calibrations/Dead Pixel Masks')
        upload_name = output_filename
        # if upload_name in sftp.listdir():
        #     print('WARNING: This will overwrite {} in pines.bu.edu:data/calibrations/Dead Pixel Masks/'.format(upload_name))
        #     upload_check = input('Do you want to continue? y/n: ')
        #     if upload_check == 'y':
        #         sftp.put(output_path,upload_name)
        #         print('Uploaded to pines.bu.edu:data/calibrations/Dead Pixel Masks/!')
        #     else:
        #         print('Skipping upload!')
        # else:
        sftp.put(output_path,upload_name)
        print('Uploaded {} to pines.bu.edu:data/calibrations/Dead Pixel Masks/!'.format(upload_name))

def dome_flat_field(date, band, lights_on_start=0, lights_on_stop=0, lights_off_start=0, lights_off_stop=0, upload=False, delete_raw=False, sftp=''):
    """Creates a master dome flat field for a given date and band. 

    :param date: date of flat images (YYYYMMDD)
    :type date: str
    :param band: band of observations (e.g., 'J' or 'H')
    :type band: str
    :param lights_on_start: start of the sequence of lights on flat images that you want to force the code to use, defaults to 0
    :type lights_on_start: int, optional
    :param lights_on_stop: end of the sequence of lights on flat images that you want to force the code to use, defaults to 0
    :type lights_on_stop: int, optional
    :param lights_off_start: start of the sequence of lights off flat images that you want to force the code to use, defaults to 0
    :type lights_off_start: int, optional
    :param lights_off_stop: end of the sequence of lights off flat images that you want to force the code to use, defaults to 0
    :type lights_off_stop: int, optional
    :param upload: whether or not to upload the resulting master flat to the PINES server, defaults to False
    :type upload: bool, optional
    :param delete_raw: whether or not to delete the raw flat files off your local machine when the flat is created, defaults to False
    :type delete_raw: bool, optional
    :param sftp: sftp connection to the PINES server, defaults to ''
    :type sftp: str, optional
    :raises RuntimeError: breaks if no suitable flats are found
    """

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
                #Download from the /data/logs/ directory on PINES.
                if not (log_path/(date+'_log.txt')).exists():
                    print('Downloading {}_log.txt to {}'.format(date,log_path))
                    sftp.get('/data/logs/'+date+'_log.txt',log_path/(date+'_log.txt'))

                log = pines_log_reader(log_path/(date+'_log.txt'))
                
                #Identify flat files. 
                flat_inds = np.where((log['Target'] == 'Flat') & (log['Filename'] != 'test.fits') & (log['Filt.'] == band))[0]
                flat_files = natsorted(list(set(log['Filename'][flat_inds]))) #Set guarantees we only grab the unique files that have been identified as flats, in case the log bugged out. 

            print('Found {} flat files.'.format(len(flat_files)))
            print('')

            #Download data to the appropriate Calibrations/Flats/Domeflats/band/Raw/ directory. 
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
        flat_files = natsorted(list(Path(dome_flat_raw_path).rglob(date+'*.fits')))
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
    lights_on_std_devs = np.zeros(num_images)
    for j in range(len(lights_on_files)):
        image_data = fits.open(dome_flat_raw_path/lights_on_files[j])[0].data[0:1024,:]
        header = fits.open(dome_flat_raw_path/lights_on_files[j])[0].header
        if header['FILTNME2'] != band:
            print('ERROR: {} taken in filter other than {}. Double check your date, try specifying start/stop file numbers, etc.'.format(lights_on_files[j], band))
            return
        flat_lights_on_cube_raw[j,:,:] = image_data #This line trims off the top two rows of the image, which are overscan.
        lights_on_std_devs[j] = np.std(image_data) #Save standard deviation of flat images to identify flats with "ski jump" horizontal bars issue.
        print(str(j+1)+'    '+str(np.mean(image_data))+'    '+str(np.std(image_data))+'    '+ str(np.std(image_data))+'    '+str(np.amin(image_data)))
    
    #Identify bad lights-on flat images (usually have bright horizontal bands at the top/bottom of images.)
    vals, lo, hi = sigmaclip(lights_on_std_devs)
    bad_locs = np.where((lights_on_std_devs < lo) | (lights_on_std_devs > hi))[0]
    good_locs = np.where((lights_on_std_devs > lo) & (lights_on_std_devs < hi))[0]
    if len(bad_locs > 0):
        print('Found {} bad lights-on flats: {}.\nExcluding from combination step.\n'.format(len(bad_locs), np.array(lights_on_files)[bad_locs]))
        time.sleep(2)
    flat_lights_on_cube_raw = flat_lights_on_cube_raw[good_locs]
    
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
    lights_off_std_devs = np.zeros(num_images)
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
        lights_off_std_devs[j] = np.std(image_data) #Save standard deviation of flat images to identify flats with "ski jump" horizontal bars issue.
        print(str(j+1)+'    '+str(np.mean(image_data))+'    '+str(np.std(image_data))+'    '+ str(np.std(image_data))+'    '+str(np.amin(image_data)))

    #Identify bad lights-on flat images (usually have bright horizontal bands at the top/bottom of images.)
    vals, lo, hi = sigmaclip(lights_off_std_devs)
    bad_locs = np.where((lights_off_std_devs < lo) | (lights_off_std_devs > hi))[0]
    good_locs = np.where((lights_off_std_devs > lo) & (lights_off_std_devs < hi))[0]
    if len(bad_locs > 0):
        print('Found {} bad lights-off flats: {}.\nExcluding from combination step.\n'.format(len(bad_locs), np.array(lights_off_files)[bad_locs]))
        time.sleep(2)
    flat_lights_off_cube_raw = flat_lights_off_cube_raw[good_locs]

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
    # if os.path.exists(output_filename):
    #     print('')
    #     print('WARNING: This will overwrite {}!'.format(output_filename))
    #     flat_check = input('Do you want to continue? y/n: ')
    #     if flat_check == 'y':
    #         hdu.writeto(output_filename,overwrite=True)
    #         print('Wrote to {}!'.format(output_filename))
    #     else:
    #         print('Not overwriting!')
    # else:
    hdu.writeto(output_filename,overwrite=True)
    print('Wrote to {}!'.format(output_filename))
    print('')

    
    if upload:
        print('Beginning upload process to pines.bu.edu...')
        print('Note, only PINES admins will be able to upload.')
        print('')
        sftp.chdir('/')
        sftp.chdir('data/calibrations/Flats/Domeflats/'+band)
        upload_name = 'master_flat_'+band+'_'+date+'.fits'
        # if upload_name in sftp.listdir():
        #     print('WARNING: This will overwrite {} in pines.bu.edu:data/calibrations/Flats/Domeflats/{}/'.format(upload_name,band))
        #     upload_check = input('Do you want to continue? y/n: ')
        #     if upload_check == 'y':
        #         sftp.put(output_filename,upload_name)
        #         print('Uploaded to pines.bu.edu:data/calibrations/Flats/Domeflats/{}/ !'.format(band))
        #     else:
        #         print('Skipping upload!')
        # else:
        sftp.put(output_filename,upload_name)
        print('Uploaded {} to pines.bu.edu:data/calibrations/Flats/Domeflats/{}/!'.format(upload_name, band))


    output_filename = pines_path/('Calibrations/Flats/Domeflats/'+band+'/Master Flats Stddev/master_flat_stddev_'+band+'_'+date+'.fits')
    #Add some header keywords detailing the master_dark creation process. 
    hdu = fits.PrimaryHDU(master_flat_error)
    hdu.header['HIERARCH DATE CREATED'] = datetime.utcnow().strftime('%Y-%m-%d')+'T'+datetime.utcnow().strftime('%H:%M:%S')

    #Now save to a file on your local machine. 
    #Check to see if other files of this name exist
    # if os.path.exists(output_filename):
    #     print('')
    #     print('WARNING: This will overwrite {}!'.format(output_filename))
    #     flat_check = input('Do you want to continue? y/n: ')
    #     if flat_check == 'y':
    #         hdu.writeto(output_filename,overwrite=True)
    #         print('Wrote to {}!'.format(output_filename))
    #     else:
    #         print('Not overwriting!')
    # else:
    hdu.writeto(output_filename,overwrite=True)
    print('Wrote to {}!'.format(output_filename))

    
    if upload:
        print('Uploading to pines.bu.edu...')
        sftp.chdir('/')
        sftp.chdir('data/calibrations/Flats Stddev/Domeflats/'+band)
        upload_name = 'master_flat_stddev_'+band+'_'+date+'.fits'
        # if upload_name in sftp.listdir():
        #     print('WARNING: This will overwrite {} in pines.bu.edu:data/calibrations/Flats/Domeflats/{}/'.format(upload_name,band))
        #     upload_check = input('Do you want to continue? y/n: ')
        #     if upload_check == 'y':
        #         sftp.put(output_filename,upload_name)
        #         print('Uploaded to pines.bu.edu:data/calibrations/Flats Stddev/Domeflats/{}/ !'.format(band))
        #     else:
        #         print('Skipping upload!')
        # else:
        sftp.put(output_filename,upload_name)
        print('Uploaded {} to pines.bu.edu:data/calibrations/Flats Stddev/Domeflats/{}/!'.format(upload_name, band))
    print('')
    if delete_raw:
        files_to_delete = glob.glob(os.path.join(dome_flat_raw_path/'*.fits'))
        for j in range(len(files_to_delete)):
            os.remove(files_to_delete[j])

    print('dome_flat_field runtime: ', np.round((time.time()-t1)/60,1), ' minutes.')
    print('Done!')

def get_calibrations(sftp, pines_path):
    """Downloads all calibration files on the PINES server. 

    :param sftp: sftp connection to the PINES server.
    :type sftp: pysftp connection
    :param pines_path: path to top-level PINES_analysis_toolkit directory containing the 'Calibrations/' folder
    :type pines_path: pathlib.PosixPath
    """

    print('Downloading calibration files from pines.bu.edu!')
    print('')
    
    sftp.chdir('/')
    sftp.chdir('data/calibrations/')
    calibration_directories = ['Bad Pixel Masks', 'Darks', 'Dead Pixel Masks', 'Flats', 'Hot Pixel Masks', 'Kokopelli Mask', 'Variable Pixel Masks']
    for i in range(len(calibration_directories)):
        sftp.chdir(calibration_directories[i])
        if calibration_directories[i] != 'Flats':
            files = sftp.listdir()
            #Set the directory the calibrations will be downloaded to.
            if calibration_directories[i] == 'Darks':
                local_path = pines_path/('Calibrations/'+calibration_directories[i]+'/Master Darks/')
            else:
                local_path = pines_path/('Calibrations/'+calibration_directories[i]+'/')

            for j in range(len(files)):
                if not (local_path/files[j]).exists():
                    sftp.get(files[j], local_path/files[j])
                    print('Downloading to {}'.format(local_path/files[j]))
        else:
            sftp.chdir('Domeflats')
            bands = sftp.listdir()
            for j in range(len(bands)):
                sftp.chdir(bands[j])
                files = sftp.listdir()
                local_path = pines_path/('Calibrations/Flats/Domeflats/'+bands[j]+'/Master Flats/')
                for k in range(len(files)):
                    if not (local_path/files[k]).exists():
                        sftp.get(files[k], local_path/files[k])
                        print('Downloading to {}'.format(local_path/files[k]))
                sftp.chdir('..')
        sftp.chdir('/')
        sftp.chdir('data/calibrations/')

def get_master_log(sftp, pines_path):
    """Downloads master_log.txt from the PINES server.

    :param sftp: sftp connection to the PINES server
    :type sftp: pysftp connection
    :param pines_path: top-level PINES path
    :type pines_path: pathlib.PosixPath
    """

    print('Downloading master_log.txt to {}/Logs/'.format(pines_path))
    sftp.chdir('/code/PINES_server_scripts/')
    sftp.get('master_log.txt',(pines_path/'Logs/master_log.txt'))
    return

def get_master_synthetic_image(sftp, target_name):
    """Downloads master synthetic image for a target from the PINES server. 

    :param sftp: sftp connection to the PINES server
    :type sftp: pysftp connection
    :param target_name: long name for the target
    :type target_name: str
    """

    pines_path = pines_dir_check()
    sftp.chdir('/data/master_images/')
    synthetic_filename = target_name.replace(' ', '')+'_master_synthetic.fits'
    download_path = pines_path/('Calibrations/Master Synthetic Images/'+synthetic_filename)
    sftp.get(synthetic_filename, download_path)
    print('Downloaded {} to {}!'.format(synthetic_filename, download_path))

def get_raw_science_files(sftp, target_name):
    """Downloads raw science files for a target

    :param sftp: sftp connection to the PINES server
    :type sftp: pysftp connection
    :param target_name: long name for the target
    :type target_name: str
    """

    t1 = time.time()
    print('')
    print('Starting get_raw_science files for {}.'.format(target_name))
   
    #Get the user's pines_analysis_toolkit path 
    pines_path = pines_dir_check()

    #Get the target's short name and set up a data directory, if necessary. 
    short_name = short_name_creator(target_name)
    if not os.path.exists(pines_path/('Objects/'+short_name)):
        object_directory_creator(pines_path, short_name)

    raw_data_path = pines_path/('Objects/'+short_name+'/raw/')
    dark_path = pines_path/('Calibrations/Darks')
    flats_path = pines_path/('Calibrations/Flats/Domeflats')
        
    #Grab an up-to-date copy of the master log, which will be used to find images. 
    get_master_log(sftp, pines_path)

    #Let's grab all of the available calibration data on pines.bu.edu.
    get_calibrations(sftp, pines_path)
    print('Calibrations up to date!')
    print('')
    time.sleep(2)

    #Read in the master target list and find images of the requested target. 
    df = pines_log_reader(pines_path/('Logs/master_log.txt'))
    targ_inds = np.where(np.array(df['Target']) == target_name)
    file_names = np.array(df['Filename'])[targ_inds]
    print('Searching pines.bu.edu for raw science files for {}.'.format(target_name))
    #Get list of dates that data are from, in chronological order. 
    dates = [int(i) for i in list(set([str.split(file_names[i],'.')[0] for i in range(len(file_names))]))]
    dates = np.array(dates)[np.argsort(dates)]
    comp_dates = np.array([int(file_names[i].split('.')[0]) for i in range(len(file_names))])
    print('Found ',len(file_names),' raw files for ',target_name,' on ',len(dates),' dates.')
    date_holder = [[] for x in range(len(dates))]
    for i in range(len(dates)):
        date = dates[i]
        print(date,': ',len(np.where(comp_dates == date)[0]),' files.')
        date_holder[i].extend(file_names[np.where(comp_dates == date)[0]])
        time.sleep(0.5)

    dates = [str(i) for i in dates]
    #Now download the identified data. 
    sftp.chdir('/data/raw/mimir')
    run_dirs = sftp.listdir()
    file_num = 1
    for i in range(len(run_dirs)):
        sftp.chdir(run_dirs[i])
        night_dirs = sftp.listdir()
        for j in range(len(night_dirs)):
            night_check = night_dirs[j]
            if night_check in dates:
                sftp.chdir(night_check)

                #Check that this night's log is in the /data/logs/ folder. If not, move a copy there.
                log_name = night_check+'_log.txt'
                if not log_name in sftp.listdir('/data/logs/'): 
                    sftp.get(log_name, pines_path/('Logs/'+log_name))
                    sftp.put(pines_path/('Logs/'+log_name), '/data/logs/'+log_name)
                
                date_holder_ind = np.where(np.array(dates) == night_check)[0][0]
                files = date_holder[date_holder_ind]
                print('Downloading {} data for {} to {}!'.format(night_check, target_name, raw_data_path))
                pbar = ProgressBar()
                for k in pbar(range(len(files))):
                    if not (raw_data_path/files[k]).exists():
                        #print('Downloading to {}, {} of {}'.format(raw_data_path/files[k],file_num,len(file_names)))
                        sftp.get(files[k],raw_data_path/files[k])
                    #else:
                        #print('{} already in {}, skipping download.'.format(files[k],raw_data_path))
                    file_num += 1
                sftp.chdir('..')
        sftp.chdir('..')
    
    print('')
    #Now grab the logs from /data/logs/ on the PINES server.
    sftp.chdir('/data/logs')
    for i in range(len(dates)):
        log_name = dates[i]+'_log.txt'
        
        print('Downloading {} to {}.'.format(log_name, pines_path/('Logs/'+log_name)))
        sftp.get(log_name, pines_path/('Logs/'+log_name))

def get_reduced_science_files(sftp, target_name):
    """Downloads reduced science files for a target from the PINES server. 

    :param sftp: sftp connection to the PINES server
    :type sftp: pysftp connection
    :param target_name: long name for the target
    :type target_name: target_name
    """

    t1 = time.time()
    
    #Get the user's pines_analysis_toolkit path 
    pines_path = pines_dir_check()

    #Get the target's short name and set up a data directory, if necessary. 
    short_name = short_name_creator(target_name)

    if not os.path.exists(pines_path/('Objects/'+short_name)):
        object_directory_creator(pines_path, short_name)

    reduced_data_path = pines_path/('Objects/'+short_name+'/reduced/')
    dark_path = pines_path/('Calibrations/Darks')
    flats_path = pines_path/('Calibrations/Flats/Domeflats')
    
    

    #Grab an up-to-date copy of the master log, which will be used to find images. 
    get_master_log(sftp, pines_path)

    #Let's grab all of the available calibration data on pines.bu.edu.
    get_calibrations(sftp, pines_path)
    print('Calibrations up to date!')
    time.sleep(2)
    
    #Read in the master target list and find images of the requested target. 
    df = pines_log_reader(pines_path/('Logs/master_log.txt'))
    targ_inds = np.where(np.array(df['Target']) == target_name)[0]
    file_names = np.array(df['Filename'])[targ_inds]
    print('')
    
    print('Searching pines.bu.edu for reduced science files for {}.'.format(target_name))
    print('')

    #Get list of dates that data are from, in chronological order. 
    dates = [int(i) for i in list(set([str.split(file_names[i],'.')[0] for i in range(len(file_names))]))]
    dates = np.array(dates)[np.argsort(dates)]
    comp_dates = np.array([int(file_names[i].split('.')[0]) for i in range(len(file_names))])
    print('Found ',len(file_names),' raw files for ',target_name,' on ',len(dates),' dates.')

    date_holder = [[] for x in range(len(dates))]
    for i in range(len(dates)):
        date = dates[i]
        print(date,': ',len(np.where(comp_dates == date)[0]),' files.')
        date_holder[i].extend(file_names[np.where(comp_dates == date)[0]])
        time.sleep(0.1)
    
    
    dates = [str(i) for i in dates]
    #Now download the identified data. 
    sftp.chdir('/data/reduced/mimir')
    run_dirs = sftp.listdir()
    file_num = 1
    for i in range(len(run_dirs)):
        sftp.chdir(run_dirs[i])
        night_dirs = sftp.listdir()
        for j in range(len(night_dirs)):
            night_check = night_dirs[j]
            if night_check in dates:
                sftp.chdir(night_check)
                date_holder_ind = np.where(np.array(dates) == night_check)[0][0]
                files = date_holder[date_holder_ind]
                files_in_path = sftp.listdir()
                for k in range(len(files)):
                    download_filename = files[k].split('.fits')[0]+'_red.fits'
                    if not (reduced_data_path/download_filename).exists():
                        if download_filename in files_in_path:
                            print('Downloading to {}, {} of {}'.format(reduced_data_path/download_filename, file_num, len(file_names)))
                            sftp.get(download_filename,reduced_data_path/download_filename)
                        else:
                            print('A reduced image does not yet exist for {}, ask an administrator to make one!'.format(files[k]))
                    else:
                        print('{} already in {}, skipping.'.format(download_filename,reduced_data_path))
                    file_num += 1
                sftp.chdir('..')
        sftp.chdir('..')

    print('')
    #Now grab the logs.
    sftp.chdir('/data/logs')
    for i in range(len(dates)):
        log_name = dates[i]+'_log.txt'
        print('Downloading {} to {}.'.format(log_name, pines_path/('Logs/'+log_name)))
        sftp.get(log_name, pines_path/('Logs/'+log_name))

    print('')
    print('get_reduced_science_files runtime: ', np.round((time.time()-t1)/60,1), ' minutes.')
    print('Done!')

def hot_pixels(date, exptime, saturation=4000., upload=False, sftp=''):
    """Creates a hot pixel mask from master dark. Flags pixels as hot if they exceed clip_lvl * the standard deviation of neighboring pixels in 
            a box of dimensions box_l x box_l surrounding the target pixel. Iterates through the master dark multiple times until no new bad pixels
            are found. 

    :param date: date the master dark was taken (YYYYMMDD)
    :type date: str
    :param exptime: exposure time of the master dark
    :type exptime: float
    :param box_l: [description], defaults to 3
    :type box_l: int, optional
    :param saturation: pixel value in ADUs over which pixels are automatically flagged as hot, defaults to 4000
    :type saturation: float, optional
    :param upload: whether or not to upload the resulting hot pixel map to the PINES server, defaults to False
    :type upload: bool, optional
    :param sftp: sftp connection to the PINES server, defaults to ''
    :type sftp: str, optional
    :raises RuntimeError: if no dark stddev files are found on disk
    """

    pines_path = pines_dir_check()
    darks_path = pines_path/('Calibrations/Darks/Master Darks/')
    all_dark_stddev_files = natsorted(list(Path(darks_path).rglob('*'+date+'.fits')))
    dark_file = ''
    for file in all_dark_stddev_files:
        if float(file.name.split('_')[2]) == exptime:
            dark_file = file
    
    if dark_file == '':
        raise RuntimeError('No dark stddev files found on disk with date '+date+' and exposure time '+str(exptime)+' seconds!')
    
    master_dark = fits.open(dark_file)[0].data
    shape = np.shape(master_dark)
    avg, med, std = sigma_clipped_stats(master_dark)

    hot_pixel_mask = np.zeros((shape[0], shape[1]), dtype='int')

    #Mask/nan any pixel that is > saturation in the master_dark
    hot_pixel_mask[np.where(master_dark > saturation)] = 1
    master_dark[np.where(master_dark > saturation)] = np.nan
    
    #Incorporate bad pixel information from the variable/Kokopelli pixel masks. 
    # variable_path = pines_path/('Calibrations/Variable Pixel Masks/vpm_'+str(exptime)+'_s_'+date+'.fits')
    # variable_mask = fits.open(variable_path)[0].data
    # master_dark[np.where(variable_mask == 1)] = np.nan
    
    # kokopelli_path = pines_path/('Calibrations/Kokopelli Mask/kokopelli_mask.fits')
    # kokopelli_mask = (1-fits.open(kokopelli_path)[0].data).astype('int')[0:1024,:]
    # master_dark[np.where(kokopelli_mask == 1)] = np.nan

    print('')
    print('Flagging hot pixels.')
    

    #This will iterate until it finds no more hot pixels. 
    total_flagged = 0
    clip_lvls = [10, 10, 10, 9, 8]
    box_ls = [13, 11, 9, 7, 5]

    iteration = 1
    for i in range(len(box_ls)):
        box_l = box_ls[i]
        clip_lvl = clip_lvls[i]
        print('Box size = {} x {}.'.format(box_l, box_l))
        print('Sigma clipping level = {}.'.format(clip_lvl))
        print('......')
        num_flagged = 999 #Initialize. Keep iterating until num_flagged = 1 or 0. 
        
        box_minus = int(box_l/2)
        box_plus = int(box_l/2) + 1
        targ_pix_ind = int(box_l**2/2)

        while (num_flagged > 1): 
            num_flagged = 0
            pbar = ProgressBar()
            for xx in pbar(range(int(box_l/2), shape[0]-int(box_l/2))):
                for yy in range(int(box_l/2), shape[1]-int(box_l/2)):
                    if (not np.isnan(master_dark[yy,xx])):
                        box_2d = master_dark[yy-box_minus:yy+box_plus,xx-box_minus:xx+box_plus] #2d cutout surrounding the target pixel.
                        box_1d = box_2d[~np.isnan(box_2d)].ravel() #Unraveled box, ignoring any NaNs. 
                        try:
                            neighbor_vals = np.delete(box_1d, targ_pix_ind) #Remove the target pixel from the 1d array. 
                            #neighbor_vals = sigmaclip(n, low=3.5, high=3.5)[0] #Clip away other bad pixels in the box to avoid biasing the mean/standard deviation.                   
                            if master_dark[yy,xx] > np.median(neighbor_vals) + clip_lvl*np.std(neighbor_vals): #Flag pixel as hot if it is more than clip_lvl sigma higher than the mean of its neighbors. 
                                hot_pixel_mask[yy,xx] = 1
                                master_dark[yy,xx] = np.nan #Set this pixel to a NaN so that it's ignored on subsequent iterations. 
                                num_flagged += 1
                        except:
                            continue

            iteration += 1
            total_flagged += num_flagged
            print('Iteration {}: {} new hot pixels identified, {} hot pixels total.'.format(iteration, num_flagged, total_flagged))
            print('')


    print('')
    print('Found {} hot pixels.'.format(total_flagged))

    output_filename = 'hpm_'+str(exptime)+'_s_'+date+'.fits'
    output_path = pines_path/('Calibrations/Hot Pixel Masks/'+output_filename)

    #Add some header keywords detailing the master_dark creation process. 
    hdu = fits.PrimaryHDU(hot_pixel_mask)
    hdu.header['HIERARCH DATE CREATED'] = datetime.utcnow().strftime('%Y-%m-%d')+'T'+datetime.utcnow().strftime('%H:%M:%S')
    hdu.header['HIERARCH SIGMA CLIP LVL'] = clip_lvl

    #Now save to a file on your local machine. 
    print('')
    print('Writing the file to '+output_filename)

    # #Check to see if other files of this name exist.
    # if os.path.exists(output_path):
    #     print('')
    #     print('WARNING: This will overwrite {}!'.format(output_path))
    #     dark_check = input('Do you want to continue? y/n: ')
    #     if dark_check == 'y':
    #         hdu.writeto(output_path,overwrite=True)
    #         print('Wrote to {}!'.format(output_path))
    #     else:
    #         print('Not overwriting!')
    # else:
    print('Wrote to {}!'.format(output_path))
    hdu.writeto(output_path,overwrite=True)
    print('')
    
    #Upload the master dark to PINES server.
    if upload:
        print('Beginning upload process to pines.bu.edu...')
        print('Note, only PINES admins will be able to upload.')
        print('')
        sftp.chdir('/')
        sftp.chdir('data/calibrations/Hot Pixel Masks')
        upload_name = output_filename
        # if upload_name in sftp.listdir():
        #     print('WARNING: This will overwrite {} in pines.bu.edu:data/calibrations/Hot Pixel Masks/'.format(upload_name))
        #     upload_check = input('Do you want to continue? y/n: ')
        #     if upload_check == 'y':
        #         sftp.put(output_path,upload_name)
        #         print('Uploaded to pines.bu.edu:data/calibrations/Hot Pixel Masks/!')
        #     else:
        #         print('Skipping upload!')
        # else:
        sftp.put(output_path,upload_name)
        print('Uploaded {} to pines.bu.edu:data/calibrations/Hot Pixel Masks/!'.format(upload_name))

def make_calibrations(sftp, exptimes, bands, dark_dates, flat_dates, dark_starts=[], dark_stops=[], lights_on_starts=[],
                      lights_on_stops=[], lights_off_starts=[], lights_off_stops=[], dark=True, flat=True, variable_pixels=True, 
                      hot_pixels=True, dead_pixels=True, bpm=True):
    """Makes all calibration files for a target. 

    :param sftp: sftp connection to the PINES server
    :type sftp: pysftp connection
    :param exptimes: list of exposure times for this target. E.g., if you observed with 15 and 30 seconds, pass [15, 30]
    :type exptimes: list of floats
    :param bands: list of bands used to observe this target. E.g., if you observed in J and H bands, pass ['J', 'H']. 
    :type bands: list of str
    :param dark_dates: list of dates on which darks were taken (order must match order of exptimes list, and there should be one entry for every entry in exptimes!)
    :type dark_dates: list of str
    :param flat_dates: list of dates on whic flats were taken (order must match order of bands list, and there should be one entry for every entry in exptimes!)
    :type flat_dates: list of str
    :param dark_starts: list of the start filenumbers for darks (order matches order of exptimes!), defaults to []
    :type dark_starts: list of int, optional
    :param dark_stops: list of the stop filenumbers for darks (order matches order of exptimes!), defaults to []
    :type dark_stops: list of int, optional
    :param lights_on_starts: list of the start filenumbers for lights on flats (order matches order of bands!), defaults to []
    :type lights_on_starts: list of int, optional
    :param lights_on_stops: list of the stop filenumbers for lights on flats (order matches order of bands!), defaults to []
    :type lights_on_stops: list of int, optional
    :param lights_off_starts: list of the start filenumbers for lights off flats (order matches order of bands!), defaults to []
    :type lights_off_starts: list of int, optional
    :param lights_off_stops: list of the stop filenumbers for lights off flats (order matches order of bands!), defaults to []
    :type lights_off_stops: list of int, optional
    :param dark: Whether or not to make the master darks, defaults to True
    :type dark: bool, optional
    :param flat: Whether or not to make the master flats, defaults to True
    :type flat: bool, optional
    :param variable_pixels: Whether or not to make the variable pixel masks, defaults to True
    :type variable_pixels: bool, optional
    :param hot_pixels: Whether or not to make the hot pixel masks, defaults to True
    :type hot_pixels: bool, optional
    :param dead_pixels: Whether or not to make the dead pixel masks, defaults to True
    :type dead_pixels: bool, optional
    :param bpm: Whether or not to make the bad pixel masks, defaults to True
    :type bpm: bool, optional
    """        

    pines_path = pat.utils.pines_dir_check()
    calibration_path = pines_path/('Calibrations/')

    assert len(bands) == len(flat_dates), 'Num. of elements in the bands list does not equal the num. of elements of the flat_dates list!\n You must have one flat_dates entry for every bands entry.'
    assert len(exptimes) == len(dark_dates), 'Num. of elements in the exptimes list does not equal the num. of elements of the dark_dates list!\n You must have one dark_dates entry for every exptimes entry.'

    #Master darks. 
    if dark:
        for i in range(len(exptimes)):
            dark_date = dark_dates[i]
            exptime = float(exptimes[i])
            dark_path = calibration_path/('Darks/Master Darks/master_dark_'+str(exptime)+'_s_'+dark_date+'.fits')
            #If the master dark for this exptime/date already exist, skip! 
            if dark_path.exists():
                print('{} already exists in {}, skipping.'.format(dark_path.name, calibration_path))
            else:
                print('Making master dark for {}-s exposure time on {}.'.format(exptime, dark_date))
                if len(dark_starts) == 0:
                    dark_start = 0
                    dark_stop  = 0
                else:
                    dark_start = dark_starts[i]
                    dark_stop  = dark_stops[i]
                pat.data.dark(dark_date, exptime, upload=True, delete_raw=True, sftp=sftp, dark_start=dark_start, dark_stop=dark_stop)

    print('')

    #Master flats. 
    if flat:
        for i in range(len(bands)):
            flat_date = flat_dates[i]
            band = bands[i]
            flat_path = calibration_path/('Flats/Domeflats/'+band+'/Master Flats/master_flat_'+band+'_'+flat_date+'.fits')
            if flat_path.exists():
                print('{} already exists in {}, skipping.'.format(flat_path.name, calibration_path))
            else:
                print('Making master flat for {}-band on {}.'.format(band, dark_date))

                if len(lights_on_starts) == 0:
                    lights_on_start  = 0
                    lights_on_stop   = 0
                    lights_off_start = 0
                    lights_off_stop  = 0
                else:
                    lights_on_start  = lights_on_starts[i]
                    lights_on_stop   = lights_on_stops[i]
                    lights_off_start = lights_off_starts[i]
                    lights_off_stop  = lights_off_stops[i]

                pat.data.dome_flat_field(flat_date, band, sftp=sftp, upload=True, delete_raw=True, lights_on_start=lights_on_start, lights_on_stop=lights_on_stop, lights_off_start=lights_off_start, lights_off_stop=lights_off_stop)

    print('')

    #Variable pixel masks. 
    if variable_pixels:
        for i in range(len(exptimes)):
            dark_date = dark_dates[i]
            exptime = float(exptimes[i])
            vpm_path = calibration_path/('Variable Pixel Masks/vpm_'+str(exptime)+'_s_'+dark_date+'.fits')
            if vpm_path.exists():
                print('{} already exists in {}, skipping.'.format(vpm_path.name, calibration_path))
            else:
                print('Making variable pixel mask for {}-s exposure time on {}.'.format(exptime, dark_date))
                pat.data.variable_pixels(dark_date, exptime, upload=True, sftp=sftp)

    print('')

    #Hot pixel masks. 
    if hot_pixels:
        for i in range(len(exptimes)):
            dark_date = dark_dates[i]
            exptime = float(exptimes[i])
            hpm_path = calibration_path/('Hot Pixel Masks/hpm_'+str(exptime)+'_s_'+dark_date+'.fits')
            if hpm_path.exists():
                print('{} already exists in {}, skipping.'.format(hpm_path.name, calibration_path))
            else:
                print('Making hot pixel mask for {}-s exposure time on {}.'.format(exptime, dark_date))
                pat.data.hot_pixels(dark_date, exptime, box_l=5, upload=True, sftp=sftp)

    print('')
    
    #Dead pixel masks
    if dead_pixels:
        for i in range(len(bands)):
            flat_date = flat_dates[i]
            band = bands[i]
            dpm_path = calibration_path/('Dead Pixel Masks/dpm_'+band+'_'+flat_date+'.fits')
            if dpm_path.exists():
                print('{} already exists in {}, skipping.'.format(dpm_path.name, calibration_path))
            else:
                print('Making dead pixel mask for {}-band exposure time on {}.'.format(band, dark_date))
                pat.data.dead_pixels(flat_date, band, upload=True, sftp=sftp)
    
    print('')
    
    #Bad pixel masks. 
    if bpm:
        for i in range(len(exptimes)):
            dark_date = dark_dates[i]
            exptime = float(exptimes[i])
            for j in range(len(bands)):
                flat_date = flat_dates[j]
                band = bands[j]
                bpm_path = calibration_path/('Bad Pixel Masks/bpm_'+band+'_'+str(exptime)+'_s_'+flat_date+'.fits')
                if bpm_path.exists():
                    print('{} already exists in {}, skipping.'.format(bpm_path.name, calibration_path))
                else:
                    print('Making dead pixel mask for {}-band, {}-s exposure time on {}.'.format(band, exptime, flat_date))
                    pat.data.bpm_maker(flat_date, dark_date, exptime, band, upload=True, sftp=sftp)
    return 

def master_dark_chooser(dark_path, header):
    """Identifies the closest dark in time to the image you want reduce, that has the correct exposure time.

    :param dark_path:  path to the ~.../PINES_analysis_toolkit/Calibrations/darks directory. 
    :type dark_path: pathlib.PosixPath
    :param header: header of the image you want to reduce
    :type header: astropy.io fits header
    :return: path to and name of the selected master dark 
    :rtype: pathlib.PosixPath, str
    """

    exptime = header['EXPTIME']
    obs_date = datetime.strptime(header['DATE-OBS'].split('T')[0].replace('-',''),'%Y%m%d')
    possible_darks = [x for x in (dark_path/'Master Darks').glob('*'+str(exptime)+'*.fits')]
    if (len(possible_darks) == 0): 
        print('ERROR: Could not find any suitable darks to reduce {}.'.format(header['FILENAME']))
        return
    else: 
        possible_dark_dates = [datetime.strptime(i.name.split('_')[-1].split('.')[0],'%Y%m%d') for i in possible_darks]
        dark_date_distances = [abs(possible_dark_dates[i]-obs_date) for i in range(len(possible_dark_dates))]
        dark_ind = np.where(np.array(dark_date_distances) == min(np.array(dark_date_distances)))[0][0]
        master_dark = fits.open(possible_darks[dark_ind])[0].data	
        master_dark_name = possible_darks[dark_ind].name
        return master_dark, master_dark_name

def master_dark_stddev_chooser(dark_std_path, header):
    """Identifies the closest master_dark_stddev in time to the image you want reduce, that has the correct exposure time.

    :param dark_std_path:  path to the ~.../PINES_analysis_toolkit/Calibrations/dark_std directory. 
    :type dark_std_path: pathlib.PosixPath
    :param header: header of the image you want to reduce
    :type header: astropy.io fits header
    :return: path to and name of the selected master dark stddev file
    :rtype: pathlib.PosixPath, str
    """

    exptime = header['EXPTIME']
    obs_date = datetime.strptime(header['DATE-OBS'].split('T')[0].replace('-',''),'%Y%m%d')
    possible_dark_stds = [x for x in dark_std_path.glob('*'+str(exptime)+'*.fits')]
    if (len(possible_dark_stds)) == 0:
        print('ERROR: Could not find any suitable master_dark_stddev images to get read noise/dark current measurement for {}.'.format(header['FILENAME']))
        return
    else:
        possible_dark_std_dates = [datetime.strptime(i.name.split('_')[-1].split('.')[0],'%Y%m%d') for i in possible_dark_stds]
        dark_std_date_distances = [abs(possible_dark_std_dates[i]-obs_date) for i in range(len(possible_dark_std_dates))]
        dark_std_ind = np.where(np.array(dark_std_date_distances) == min(np.array(dark_std_date_distances)))[0][0]
        master_dark_std = fits.open(possible_dark_stds[dark_std_ind])[0].data
        master_dark_std_name = possible_dark_stds[dark_std_ind].name
        return master_dark_std

def master_flat_chooser(flats_path, header):
    """Identifies the closest master flat in time to the image you want reduce, that has the correct exposure time.

    :param flats_path:  path to the ~.../PINES_analysis_toolkit/Calibrations/Flats/Domeflats/ directory. 
    :type dark_path: pathlib.PosixPath
    :param header: header of the image you want to reduce
    :type header: astropy.io fits header
    :return: path to and name of the selected master flat
    :rtype: pathlib.PosixPath, str
    """
    band = header['FILTNME2']
    obs_date = datetime.strptime(header['DATE-OBS'].split('T')[0].replace('-',''),'%Y%m%d')
    possible_flats = [x for x in (flats_path/(band+'/Master Flats')).glob('*.fits')]
    if (len(possible_flats)) == 0:
        print('ERROR: Could not find any suitable flats to reduce {}.'.format(header['FILENAME']))
        return
    else:
        possible_flat_dates = [datetime.strptime(i.name.split('_')[-1].split('.')[0],'%Y%m%d') for i in possible_flats]
        flat_date_distances = [abs(possible_flat_dates[i]-obs_date) for i in range(len(possible_flat_dates))]
        flat_ind = np.where(np.array(flat_date_distances) == min(np.array(flat_date_distances)))[0][0]
        master_flat = fits.open(possible_flats[flat_ind])[0].data
        master_flat_name = possible_flats[flat_ind].name
        return master_flat, master_flat_name

def reduce(short_name, upload=False, delete_raw=False, delete_reduced=False, sftp='', manual_flat_path='', manual_dark_path='', manual_bpm_path='', linearity_correction=False, force_output_path=''):
	"""Reduces raw PINES science images and writes them out to disk.

	:param short_name: the short name for the target
	:type short_name: str
	:param upload:  whether or not to upload reduced images to the PINES server, defaults to False
	:type upload: bool, optional
	:param delete_raw: whether or not to delete local raw images when done making reduced images, defaults to False
	:type delete_raw: bool, optional
	:param delete_reduced: whether or not to delete local reduced images after upload to the server, defaults to False
	:type delete_reduced: bool, optional
	:param sftp: sftp connection to the PINES server, defaults to ''
	:type sftp: str, optional
	:param manual_flat_path: path to dome flat you want to force the reduction to use, defaults to ''
	:type manual_flat_path: pathlib.PosixPath, optional
	:param manual_dark_path: path to dark you want to force the reduction to use, defaults to ''
	:type manual_dark_path: pathlib.PosixPath, optional
	:param manual_bpm_path: path to bad pixel mask want to force the reduction to use, defaults to ''
	:type manual_bpm_path: pathlib.PosixPath, optional
	:param linearity_correction: whether or not to try a linearity correction, defaults to False
	:type linearity_correction: bool, optional
	:param force_output_path: user-chosen directory to use in place of the default ~/Documents/PINES_analysis_toolkit directory for analysis, defaults to ''
	:type force_output_path: str, optional
	"""

	t1 = time.time()
	print('')
	if (upload is True) and (sftp == ''):
		print('ERROR: You must pass an sftp connection if you want to upload reduced files to pines.bu.edu.!')
		return

	if force_output_path != '':
		pines_path = force_output_path
	else:
		pines_path = pines_dir_check()
		

	#Paths
	raw_path = pines_path/('Objects/'+short_name+'/raw')
	raw_files = [Path(i) for i in natsorted(glob.glob(os.path.join(raw_path,'*.fits')))] #Natsort sorts things how you expect them to be sorted.
	dark_path = pines_path/('Calibrations/Darks/')
	reduced_path = pines_path/('Objects/'+short_name+'/reduced')
	flats_path = pines_path/('Calibrations/Flats/Domeflats')
	bpm_path = pines_path/('Calibrations/Bad Pixel Masks')

	#Now begin the loop to load and reduce the raw science data
	print("Reducing data for {}.".format(short_name))
	print('Note: any reduced data already in reduced directory will be not be re-reduced.')
				
	pbar = ProgressBar()
	for i in pbar(range(np.size(raw_files))):
		target_filename = reduced_path/(raw_files[i].name.split('.fits')[0]+'_red.fits')
		if target_filename.exists():
			#print('{} already in reduced directory, skipping.'.format(target_filename.name))
			continue
		hdulist = fits.open(raw_files[i])
		header = hdulist[0].header

		try:
			frame_raw = fits.open(raw_files[i])[0].data.astype('float32')
		except:
			pdb.set_trace()


		frame_raw = frame_raw[0:1024,:] #Cuts off 2 rows of overscan (?) pixels
			
		#Add a flag to the header to check if background is near saturation
		if sigma_clipped_stats(frame_raw)[1] > 3800: 
			sat_flag = 1
		else:
			sat_flag = 0

		#Load in the dark/flat files. If a manual path is not provided, choose the reduction image that is closest in time to when the image was taken.
		if manual_flat_path == '':
			master_flat, master_flat_name = master_flat_chooser(flats_path, header)
		else:
			master_flat = fits.open(manual_flat_path)[0].data
			master_flat_name = manual_flat_path.name

		if manual_dark_path == '':
			master_dark, master_dark_name = master_dark_chooser(dark_path, header)
		else:
			master_dark = fits.open(manual_dark_path)[0].data
			master_dark_name = manual_dark_path.name
		
		if manual_bpm_path == '':
			bad_pixel_mask, bad_pixel_mask_name = bpm_chooser(bpm_path, header)
		else:
			bad_pixel_mask = fits.open(manual_bpm_path)[0].data
			bad_pixel_mask_name = manual_bpm_path.name

		#If linearity_correction is set, correct the pixels for linearity.
		if linearity_correction:
			log_line_coeffs_path = pines_path/('Calibrations/Linearity/log_line_coeffs.fits')
			hdu_list = fits.open(log_line_coeffs_path)
			log_line_coeffs = hdu_list[0].data[:,0:1024,:]
			median_linear_count_rate = np.nanmedian(10**log_line_coeffs[0,:,:])	
			frame_raw[np.where(bad_pixel_mask == 1)] = np.nan
			measured_counts = frame_raw - master_dark
			corrected_counts = (measured_counts/(10**log_line_coeffs[0,:,:]))**(1/log_line_coeffs[1,:,:]) * median_linear_count_rate
			frame_red = corrected_counts
			breakpoint()
		else:
			#Reduce the image. 
			frame_red = (frame_raw - master_dark)/master_flat
			frame_red = frame_red.astype('float32')
			
			#Set bad pixels to NaNs. 
			frame_red[np.where(bad_pixel_mask == 1)] = np.nan
		
		# #Do a background model subtraction
		# frame_red = bg_2d(frame_red)

		#Store some parameters in the reduced file's header. 
		#Naming convention follows https://docs.astropy.org/en/stable/io/fits/usage/headers.html.
		header['HIERARCH DATE REDUCED'] = datetime.utcnow().strftime('%Y-%m-%d')+'T'+datetime.utcnow().strftime('%H:%M:%S')
		header['HIERARCH MASTER DARK'] = master_dark_name
		header['HIERARCH MASTER FLAT'] = master_flat_name
		header['HIERARCH BAD PIXEL MASK'] = bad_pixel_mask_name
		header['HIERARCH SATURATION FLAG'] = sat_flag

		if not os.path.exists(target_filename):
			fits.writeto(target_filename, frame_red, header)
			#print("Reducing {}: {} of {}, band = {}, exptime = {} s, dark = {}, flat = {}".format(raw_files[i].name, str(i+1), str(np.size(raw_files)), header['FILTNME2'], header['EXPTIME'], master_dark_name, master_flat_name))
				
	if upload:
		print('')
		print('Beginning upload process to pines.bu.edu...')
		print('NOTE:    Only PINES admins are able to upload.')
		print('WARNING: If these reduced images already exist on the PINES server, they will be overwritten!')
		time.sleep(1)
		sftp.chdir('/data/reduced/mimir/')
		files_to_upload = np.array(natsorted(np.array([x for x in reduced_path.glob('*.fits')])))
		print('Uploading reduced {} data to the PINES server!'.format(short_name))
		pbar = ProgressBar()
		for i in pbar(range(len(files_to_upload))):
			file = files_to_upload[i]
			night_name = files_to_upload[i].name.split('.')[0]
			for dir_change in sftp.listdir():
				sftp.chdir(dir_change)
				nights = sftp.listdir()
				if night_name in nights:
					ind  = np.where(np.array(nights) == night_name)[0][0]
					break
				sftp.chdir('..')
			if nights[ind] != night_name:
				print('ERROR: the date of the file you want to upload does not match the date directory where the program wants to upload it.')
				pdb.set_trace()
			
			#print('Uploading to {}/{}, {} of {}'.format(sftp.getcwd(),nights[ind]+'/'+file.name, i+1,len(files_to_upload)))
			sftp.put(file,nights[ind]+'/'+file.name)
			sftp.chdir('..')
			

	if delete_raw:
		files_to_delete = glob.glob(os.path.join(raw_path/'*.fits'))
		for j in range(len(files_to_delete)):
			os.remove(files_to_delete[j])

	if delete_reduced:
		files_to_delete = glob.glob(os.path.join(reduced_path/'*.fits'))
		for j in range(len(files_to_delete)):
			os.remove(files_to_delete[j])

def straggler_upload(sftp, targets):
    """Uploads images of targets that were not previously logged.

    :param sftp: sftp connection to the PINES server
    :type sftp: pysftp connection
    :param targets: list of target full names you want to look for stragglers to upload
    :type targets: list of str
    """

    sftp.chdir('/data/reduced/mimir/')
    runs = sftp.listdir()
    pines_path = pines_dir_check()
    for i in range(len(targets)):
        target = targets[i]
        short_name = short_name_creator(target)
        raw_path = pines_path/('Objects/'+short_name+'/raw/')
        red_path = pines_path/('Objects/'+short_name+'/reduced/')
        raw_files = np.array(natsorted(glob(str(raw_path)+'/*.fits')))
        for j in range(len(raw_files)):
            sftp.chdir('/data/reduced/mimir/')
            red_file = str(red_path)+'/'+raw_files[j].split('/')[-1].split('.fits')[0]+'_red.fits'
            night = red_file.split('/')[-1].split('.')[0]
            run_guess = night[0:6]
            ind = np.where(np.array(runs) == run_guess)[0][0]
            if ind + 1 == len(runs):
                inds = np.arange(ind-1, ind+1)
                runs = np.array(runs)[inds]
                runs = [runs[1], runs[0]]
            else:
                inds = np.arange(ind-1, ind+2)
                runs = np.array(runs)[inds]
                runs = [runs[1], runs[0], runs[2]] 
            
            found = False
            for k in range(len(runs)):
                sftp.chdir(runs[k])
                if night in sftp.listdir():
                    print('Uploading {} to {}.'.format(red_file.split('/')[-1], sftp.pwd+'/'+night))
                    sftp.chdir(night)
                    sftp.put(red_file, red_file.split('/')[-1])
                    found = True
                else:
                    sftp.chdir('..')
                if found:
                    break

def upload_reduced_data(sftp, short_name): 
    """Standalone function for uploading reduced data to the PINES server

    :param sftp: sftp connection to the PINES server
    :type sftp: pysftp connection
    :param short_name: short name for the target
    :type short_name: str 
    """
    pines_path = pat.utils.pines_dir_check()
    reduced_path = pines_path/('Objects/'+short_name+'/reduced')

    print('')
    print('Beginning upload process to pines.bu.edu...')
    print('NOTE:    Only PINES admins are able to upload.')
    print('WARNING: If these reduced images already exist on the PINES server, they will be overwritten!')

    sftp.chdir('/data/reduced/mimir/')
    files_to_upload = np.array(natsorted(np.array([x for x in reduced_path.glob('*red.fits')])))
    print('Uploading reduced {} data to the PINES server!'.format(short_name))
    pbar = ProgressBar()
    for i in pbar(range(len(files_to_upload))):
        file = files_to_upload[i]
        night_name = files_to_upload[i].name.split('.')[0]
        for dir_change in sftp.listdir():
            sftp.chdir(dir_change)
            nights = sftp.listdir()
            if night_name in nights:
                ind  = np.where(np.array(nights) == night_name)[0][0]
                break
            sftp.chdir('..')
        if nights[ind] != night_name:
            print('ERROR: the date of the file you want to upload does not match the date directory where the program wants to upload it.')
            pdb.set_trace()
        
        #print('Uploading to {}/{}, {} of {}'.format(sftp.getcwd(),nights[ind]+'/'+file.name, i+1,len(files_to_upload)))
        sftp.put(file,nights[ind]+'/'+file.name)
        sftp.chdir('..')
        
def variable_pixels(date, exptime, clip_lvl=5., upload=False, sftp=''):
    """Creates a map of variable pixels using a dark stddev image

    :param date: UT date during which the dome flat field data was obtained (YYYYMMDD)
    :type date: str
    :param exptime: exposure time of the dark images in question
    :type exptime: float
    :param clip_lvl: sigma level above which pixels are considered to be variable, defaults to 5
    :type clip_lvl: float, optional
    :param upload: whether or not to upload the variable pixel map to the PINES server, defaults to False
    :type upload: bool, optional
    :param sftp: sftp connection to the PINES server, defaults to ''
    :type sftp: pysftp connection, optional
    :raises RuntimeError: if no dark stddev files are found on disk
    """
    pines_path = pines_dir_check()
    dark_stddev_path = pines_path/('Calibrations/Darks/Master Darks Stddev/')
    all_dark_stddev_files = natsorted(list(Path(dark_stddev_path).rglob('*'+date+'.fits')))
    dark_stddev_file = ''
    for file in all_dark_stddev_files:
        if float(file.name.split('_')[3]) == exptime:
            dark_stddev_file = file

    if dark_stddev_file == '':
        raise RuntimeError('No dark stddev files found on disk with date '+date+' and exposure time '+str(exptime)+' seconds!')
    
    master_dark_stddev = fits.open(dark_stddev_file)[0].data
    shape = np.shape(master_dark_stddev)
    #Flag variable pixels, those that vary >= clip_lvl * the mean variation. Save a boolean mask of variable pixels to incorporate into the bad pixel mask. 

    variable_inds = np.where(master_dark_stddev >= clip_lvl*np.nanmean(master_dark_stddev))
    variable_mask = np.zeros((shape[0], shape[1]), dtype='int')
    variable_mask[variable_inds] = 1
    print('')
    print('Found {} variable pixels.'.format(len(variable_inds[0])))

    output_filename = 'vpm_'+str(exptime)+'_s_'+date+'.fits'
    output_path = pines_path/('Calibrations/Variable Pixel Masks/'+output_filename)

    #Add some header keywords detailing the master_dark creation process. 
    hdu = fits.PrimaryHDU(variable_mask)
    hdu.header['HIERARCH DATE CREATED'] = datetime.utcnow().strftime('%Y-%m-%d')+'T'+datetime.utcnow().strftime('%H:%M:%S')
    hdu.header['HIERARCH SIGMA CLIP LVL'] = clip_lvl
    hdu.header['HIERARCH CLIP VALUE'] = clip_lvl*np.nanmean(master_dark_stddev)

    #Now save to a file on your local machine. 
    print('')
    # #Check to see if other files of this name exist.
    # if os.path.exists(output_path):
    #     print('')
    #     print('WARNING: This will overwrite {}!'.format(output_path))
    #     dark_check = input('Do you want to continue? y/n: ')
    #     if dark_check == 'y':
    #         hdu.writeto(output_path,overwrite=True)
    #         print('Wrote to {}!'.format(output_path))
    #     else:
    #         print('Not overwriting!')
    # else:
    hdu.writeto(output_path,overwrite=True)
    print('Wrote to {}!'.format(output_path))

    #Upload the master dark to PINES server.
    if upload:
        print('Uploading to pines.bu.edu...')
        sftp.chdir('/')
        sftp.chdir('data/calibrations/Variable Pixel Masks')
        upload_name = output_filename
        # if upload_name in sftp.listdir():
        #     print('WARNING: This will overwrite {} in pines.bu.edu:data/calibrations/Variable Pixel Masks/'.format(upload_name))
        #     upload_check = input('Do you want to continue? y/n: ')
        #     if upload_check == 'y':
        #         sftp.put(output_path,upload_name)
        #         print('Uploaded to pines.bu.edu:data/calibrations/Variable Pixel Masks/!')
        #     else:
        #         print('Skipping upload!')
        # else:
        sftp.put(output_path,upload_name)
        print('Uploaded {} to pines.bu.edu:data/calibrations/Variable Pixel Masks/!'.format(upload_name))

