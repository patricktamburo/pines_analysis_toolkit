import pdb
from pathlib import Path
import numpy, scipy
from scipy.integrate import *
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
import subprocess 
import os, glob
import sys, math
from scipy import stats, signal
from datetime import datetime
from numpy import *
from scipy import optimize
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
import time
import multiprocessing
from multiprocessing import Pool
import pickle
from astropy.io import fits
import pdb
from pathlib import Path
from astropy.stats import sigma_clipped_stats
import matplotlib.pyplot as plt
import natsort
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.data.get_master_dome_flats import get_master_dome_flats
from pines_analysis_toolkit.data.get_master_darks import get_master_darks
import getpass
import paramiko

'''Authors:
		Paul Dalba, Boston University, February 2017
		Patrick Tamburo, Boston University, July 2019 and June 2020
	Purpose:
		Reduces raw PINES science images for a specified target and saves the reduced images. 
	Inputs:
		target_name (str): the target's full 2MASS name, i.e. '2MASS J01234567+0123456' 
		flat_type (str, optional): either 'dome' or 'sky'. By default, set to dome. We haven't seen better performance with sky flats. 
		upload (bool, optional): whether or not to upload the reduced images to pines.bu.edu. By default, False (so you won't try to upload!).
	Outputs:
		None
	TODO:
		Change things to Path objects...current approach of splitting path names on '/' will break on Windows.
		Protection for overwriting files?
'''
	

def reduce(target_name, flat_type='dome', upload=False):

	def prev_month(run_string):
		#If a calibration file is not found with the current run_string, change run_string to be one month before
		#and try to restore those calibrations. 
		year = int(run_string[0:4])
		month = int(run_string[4:])
		if month != 1:
			month = month - 1
			if month < 10:
				new_month_string = '0'+str(month)
			else:
				new_month_string = str(month)
			new_run_string = str(year)+new_month_string
		else:
			new_month_string = '12'
			year = year - 1
			new_year_string = str(year)
			new_run_string = new_year_string + new_month_string

		return new_run_string
	
	t1 = time.time()

	pines_path = pines_dir_check()
	short_name = short_name_creator(target_name)

	#Paths
	#TODO fix these!
	raw_path = pines_path+'Objects/'+short_name+'/raw/'
	raw_files = natsort.natsorted(glob.glob(raw_path+'*.fits')) #Natsort sorts things how you expect them to be sorted.
	dark_path = pines_path+'Calibrations/Darks/'
	flats_path = pines_path+'Calibrations/Flats/'
	if flat_type == 'sky':
		reduced_path = pines_path+'Objects/'+short_name+'/skyflat_reduced/'
	elif flat_type == 'dome':
		reduced_path = pines_path+'Objects/'+short_name+'/reduced/'
		
	#Prompt login: 
	username = input('Enter username: ')
	password = getpass.getpass('Enter password: ')

	t1 = time.time()

	#Open ssh connection and set up local/remote paths.
	ssh = paramiko.SSHClient()
	ssh.load_system_host_keys()
	ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
	ssh.connect('pines.bu.edu',username=username, password=password)
	password = ''
	sftp = ssh.open_sftp()
	
	#Let's grab all of the available calibration data on pines.bu.edu.
	print('')
	get_master_dome_flats(sftp, pines_path+'Calibrations/Flats/Domeflats/')
	get_master_darks(sftp, pines_path+'Calibrations/Darks/')
	print('Domeflats and darks up to date!')
	print('')
	time.sleep(2)

	#Now begin the loop to load and reduce the raw science data
	for i in range(size(raw_files)):
		target_filename = reduced_path+'.'.join(raw_files[i].split('/')[-1].split('.')[:-1])+'_red.fits'
		if os.path.isfile(target_filename):
			print(raw_files[i].split('/')[-1]+' already in reduced directory, skipping.')
			continue
		hdulist = fits.open(raw_files[i])
		header = hdulist[0].header
		
		frame_raw = fits.open(raw_files[i])[0].data

		frame_raw = frame_raw[0:1024,:] #Cuts off 2 rows of overscan (?) pixels

		#Get the band and exptime to load in the appropriate dark/flat files.
		band = header['FILTNME2']
		exptime = header['EXPTIME']
		obs_date = datetime.strptime(header['DATE-OBS'].split('T')[0].replace('-',''),'%Y%m%d')

		possible_darks = glob.glob(pines_path+'Calibrations/Darks/Master Darks/*['+str(exptime)+']*.fits')
		possible_flats = glob.glob(pines_path+'Calibrations/Flats/Domeflats/'+band+'/Master Flats/*.fits')
		
		if (len(possible_darks) == 0): 
			print('ERROR: Could not find any suitable darks to reduce {}. Check exposure time of the image, and how the possible_darks list is made.'.format(raw_files[i].split('/')[-1]))
			pdb.set_trace()
		if (len(possible_flats)) == 0:
			print('ERROR: Could not find any suitable flatsto reduce {}. Check filter of the image, and how the possible_flats list is made.'.format(raw_files[i].split('/')[-1]))
			pdb.set_trace()
		
		possible_dark_dates = [datetime.strptime(possible_darks[i].split('_')[-1].split('.')[0],'%Y%m%d') for i in range(len(possible_darks))]
		possible_flat_dates = [datetime.strptime(possible_flats[i].split('_')[-1].split('.')[0],'%Y%m%d') for i in range(len(possible_flats))]

		dark_date_distances = [abs(possible_dark_dates[i]-obs_date) for i in range(len(possible_dark_dates))]
		dark_ind = where(array(dark_date_distances) == min(array(dark_date_distances)))[0][0]
		master_dark = fits.open(possible_darks[dark_ind])[0].data	
		master_dark_name = possible_darks[dark_ind].split('/')[-1]

		if flat_type == 'dome':  
			flat_date_distances = [abs(possible_flat_dates[i]-obs_date) for i in range(len(possible_flat_dates))]
			flat_ind = where(array(flat_date_distances) == min(array(flat_date_distances)))[0][0]
			master_flat = fits.open(possible_flats[flat_ind])[0].data
			master_flat_name = possible_flats[flat_ind].split('/')[-1]

		elif flat_type == 'sky':
			print('ERROR: Still have to add skyflat functionality.')
			pdb.set_trace()
		
		frame_red = (frame_raw - master_dark)/master_flat

		#Store some parameters in the reduced file's header. 
		#Naming convention follows https://docs.astropy.org/en/stable/io/fits/usage/headers.html.
		header['HIERARCH REDUCER'] = (username, 'The pines.bu.edu username of the user who reduced this image')
		header['HIERARCH DATE REDUCED'] = (datetime.utcnow().strftime('%Y-%m-%d')+'T'+datetime.utcnow().strftime('%H:%M:%S'), 'The UTC date on which the image was reduced')
		header['HIERARCH MASTER DARK'] = (master_dark_name, 'The name of the master_dark image used to do the reduction')
		header['HIERARCH MASTER FLAT'] = (master_flat_name, 'The name of the master_flat image used to do the reduction')
		header['HIERARCH FLAT TYPE'] = (flat_type, "'dome' or 'sky'")
		if not os.path.exists(target_filename):
			fits.writeto(target_filename, frame_red, header)
			print('')
			print("Reducing {}: {} of {}, band = {}, exptime = {} s, dark = '{}', flat = '{}'".format(raw_files[i].split('/')[-1], str(i+1), str(size(raw_files)), band, str(exptime), master_dark_name, master_flat_name))
		else:
			print('{} already in reduced path, skipping.'.format(raw_files[i].split('/')[-1]))
	
	print('')
	username = ''
	if upload:
		print('Beginning upload process to pines.bu.edu...')
		print('Note, only PINES admins will be able to upload.')
		print('')
		time.sleep(2)
		sftp.chdir('data/reduced/mimir/')
		files_to_upload = natsort.natsorted(glob.glob(reduced_path+'*_red.fits'))
		for i in range(len(files_to_upload)):
			file = files_to_upload[i]
			night_name = file.split('/')[-1].split('.')[0]

			for dir_change in sftp.listdir():
				sftp.chdir(dir_change)
				nights = sftp.listdir()
				if night_name in nights:
					ind  = where(array(nights) == night_name)[0][0]
					break
				sftp.chdir('..')
			if nights[ind] != night_name:
				print('ERROR: the date of the file you want to upload does not match the date directory where the program wants to upload it.')
				pdb.set_trace()
			if file.split('/')[-1] not in sftp.listdir(nights[ind]):
				print('Uploading to {}/{}, {} of {}'.format(sftp.getcwd(),nights[ind]+'/'+file.split('/')[-1], i+1,len(files_to_upload)))
				sftp.put(file,nights[ind]+'/'+file.split('/')[-1])
				sftp.chdir('..')
			else:
				print('{} already exists in {}, skipping.'.format(file.split('/')[-1],nights[ind]))
				sftp.chdir('..')

		
	print('reduce runtime: ', np.round((time.time()-t1)/60,1), ' minutes.')
	print('Reduction script has completed.')
	print('')