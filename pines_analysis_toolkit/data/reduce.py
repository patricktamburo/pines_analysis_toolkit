import pdb
from pathlib import Path
import scipy
from scipy.integrate import *
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
import subprocess 
import os, glob
import sys, math
from scipy import stats, signal
from datetime import datetime
import numpy as np
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
from pines_analysis_toolkit.data.master_flat_chooser import master_flat_chooser
from pines_analysis_toolkit.data.master_dark_chooser import master_dark_chooser
from pines_analysis_toolkit.data.bpm_chooser import bpm_chooser
import getpass
import paramiko
import pysftp 
from pines_analysis_toolkit.utils.quick_plot import quick_plot as qp
from pines_analysis_toolkit.data.bg_2d import bg_2d

'''Authors:
		Paul Dalba, Boston University, February 2017
		Patrick Tamburo, Boston University, July 2019 and June 2020
	Purpose:
		Reduces raw PINES science images for a specified target and saves the reduced images. 
	Inputs:
		target_name (str): the target's full 2MASS name, i.e. '2MASS J01234567+0123456' 
		upload (bool, optional): whether or not to upload the reduced images to pines.bu.edu. By default, False (so you won't try to upload!).
		sftp (pysftp.Connection, optional): the sftp connection to the pines server, required if you are going to upload reduced data
		delete_raw (bool, optional): whether or not to delete raw files from your local directory for this target when done reduction/upload process. 
		delete_reduced (bool, optional): whether or not to delete reduced files from your local directory for this target when done reduction/upload process.
	Outputs:
		None
	TODO:
		None
'''
	
def reduce(target_name, upload=False, delete_raw=False, delete_reduced=False, sftp=''):
	t1 = time.time()
	print('')
	print('Starting reduction script for {}.'.format(target_name))

	if (upload is True) and (sftp == ''):
		print('ERROR: You must pass an sftp connection if you want to upload reduced files to pines.bu.edu.!')
		return

	pines_path = pines_dir_check()
	short_name = short_name_creator(target_name)

	#Paths
	raw_path = pines_path/('Objects/'+short_name+'/raw')
	raw_files = [Path(i) for i in natsort.natsorted(glob.glob(os.path.join(raw_path,'*.fits')))] #Natsort sorts things how you expect them to be sorted.
	dark_path = pines_path/('Calibrations/Darks')
	reduced_path = pines_path/('Objects/'+short_name+'/reduced')
	flats_path = pines_path/('Calibrations/Flats/Domeflats')
	bpm_path = pines_path/('Calibrations/Bad Pixel Masks')

	#Now begin the loop to load and reduce the raw science data
	for i in range(np.size(raw_files)):
		target_filename = reduced_path/(raw_files[i].name.split('.fits')[0]+'_red.fits')
		if target_filename.exists():
			print('{} already in reduced directory, skipping.'.format(target_filename.name))
			continue
		hdulist = fits.open(raw_files[i])
		header = hdulist[0].header
		
		frame_raw = fits.open(raw_files[i])[0].data

		frame_raw = frame_raw[0:1024,:] #Cuts off 2 rows of overscan (?) pixels

		#Load in the dark/flat files that were taken as close in time as possible to the target image. 
		master_flat, master_flat_name = master_flat_chooser(flats_path, header)
		master_dark, master_dark_name = master_dark_chooser(dark_path, header)
		bad_pixel_mask, bad_pixel_mask_name = bpm_chooser(bpm_path, header)

		#Reduce the image. 
		frame_red = (frame_raw - master_dark)/master_flat
		frame_red = frame_red.astype('float32')

		#Set bad pixels to NaNs. 
		frame_red[np.where(bad_pixel_mask == 1)] = np.nan
		#qp(frame_red)

		# #Do a background model subtraction
		# frame_red = bg_2d(frame_red)

		#Store some parameters in the reduced file's header. 
		#Naming convention follows https://docs.astropy.org/en/stable/io/fits/usage/headers.html.
		header['HIERARCH DATE REDUCED'] = datetime.utcnow().strftime('%Y-%m-%d')+'T'+datetime.utcnow().strftime('%H:%M:%S')
		header['HIERARCH MASTER DARK'] = master_dark_name
		header['HIERARCH MASTER FLAT'] = master_flat_name
		header['HIERARCH BAD PIXEL MASK'] = bad_pixel_mask_name

		if not os.path.exists(target_filename):
			fits.writeto(target_filename, frame_red, header)
			print('')
			print("Reducing {}: {} of {}, band = {}, exptime = {} s, dark = {}, flat = {}".format(raw_files[i].name, str(i+1), str(np.size(raw_files)), header['FILTNME2'], header['EXPTIME'], master_dark_name, master_flat_name))
		else:
			print('{} already in reduced path, skipping.'.format(raw_files[i].name))
	
	print('')
	if upload:
		print('Beginning upload process to pines.bu.edu...')
		print('Note, only PINES admins will be able to upload.')
		print('')
		time.sleep(2)
		sftp.chdir('/data/reduced/mimir/')
		files_to_upload = array(natsort.natsorted(array([x for x in reduced_path.glob('*.fits')])))
		for i in range(len(files_to_upload)):
			file = files_to_upload[i]
			night_name = files_to_upload[i].name.split('.')[0]
			for dir_change in sftp.listdir():
				sftp.chdir(dir_change)
				nights = sftp.listdir()
				if night_name in nights:
					ind  = np.where(array(nights) == night_name)[0][0]
					break
				sftp.chdir('..')
			if nights[ind] != night_name:
				print('ERROR: the date of the file you want to upload does not match the date directory where the program wants to upload it.')
				pdb.set_trace()
			if file.name not in sftp.listdir(nights[ind]):
				print('Uploading to {}/{}, {} of {}'.format(sftp.getcwd(),nights[ind]+'/'+file.name, i+1,len(files_to_upload)))
				sftp.put(file,nights[ind]+'/'+file.name)
				sftp.chdir('..')
			else:
				print('{} already exists in {}, skipping. {} of {}.'.format(file.name,nights[ind],i+1,len(files_to_upload)))
				sftp.chdir('..')

	if delete_raw:
		files_to_delete = glob.glob(os.path.join(raw_path/'*.fits'))
		for j in range(len(files_to_delete)):
			os.remove(files_to_delete[j])

	if delete_reduced:
		files_to_delete = glob.glob(os.path.join(reduced_path/'*.fits'))
		for j in range(len(files_to_delete)):
			os.remove(files_to_delete[j])

	print('reduce runtime: ', round((time.time()-t1)/60,1), ' minutes.')
	print('Reduction script has completed.')
	print('')