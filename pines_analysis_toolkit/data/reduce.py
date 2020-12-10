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
from progressbar import ProgressBar
from astropy.visualization import ImageNormalize, ZScaleInterval, SquaredStretch, SqrtStretch, SinhStretch
from copy import deepcopy
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
		Add check to see if .fits files are the correct size. If they didn't fully download, fits.open() will crash. 
'''
	
def reduce(target_name, upload=False, delete_raw=False, delete_reduced=False, sftp=''):
	t1 = time.time()
	print('')
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
	print("Reducing data for {}.".format(target_name))
	print('Note: any reduced data already in reduced directory will be not be re-reduced.')
				
	pbar = ProgressBar()
	for i in pbar(range(np.size(raw_files))):
		target_filename = reduced_path/(raw_files[i].name.split('.fits')[0]+'_red.fits')
		if target_filename.exists():
			#print('{} already in reduced directory, skipping.'.format(target_filename.name))
			continue
		hdulist = fits.open(raw_files[i])
		header = hdulist[0].header

		frame_raw = fits.open(raw_files[i])[0].data


		frame_raw = frame_raw[0:1024,:] #Cuts off 2 rows of overscan (?) pixels

		#Add a flag to the header to check if background is near saturation
		if sigma_clipped_stats(frame_raw)[1] > 3800: 
			sat_flag = 1
		else:
			sat_flag = 0

		#Load in the dark/flat files that were taken as close in time as possible to the target image. 
		master_flat, master_flat_name = master_flat_chooser(flats_path, header)
		master_dark, master_dark_name = master_dark_chooser(dark_path, header)
		bad_pixel_mask, bad_pixel_mask_name = bpm_chooser(bpm_path, header)

		#Reduce the image. 
		frame_red = (frame_raw - master_dark)/master_flat
		frame_red = frame_red.astype('float32')
		
		#frame_red_no_flag = deepcopy(frame_red)

		#Set bad pixels to NaNs. 
		frame_red[np.where(bad_pixel_mask == 1)] = np.nan

		# norm = ImageNormalize(data=frame_red, interval=ZScaleInterval())
		# plt.ion()
		# fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(16,7), sharex=True, sharey=True)
		# im = ax[0].imshow(frame_red_no_flag, origin='lower', norm=norm)
		# ax[0].set_title('Reduced')
		# divider = make_axes_locatable(ax[0])
		# cax = divider.new_vertical(size="3%", pad=0.3, pack_start=True)
		# fig.add_axes(cax)
		# fig.colorbar(im, cax=cax, orientation="horizontal", label='ADUs')

		# im2 = ax[1].imshow(frame_red, origin='lower', norm=norm)
		# ax[1].set_title('Reduced, Bad Pixels Flagged')
		# divider = make_axes_locatable(ax[1])
		# cax = divider.new_vertical(size="3%", pad=0.3, pack_start=True)
		# fig.add_axes(cax)
		# fig.colorbar(im2, cax=cax, orientation="horizontal", label='ADUs')

		# kernel = Gaussian2DKernel(x_stddev=0.5)
		# frame_red_interp = interpolate_replace_nans(frame_red, kernel=kernel)
		# im3 = ax[2].imshow(frame_red_interp, origin='lower', norm=norm)
		# ax[2].set_title('Reduced, Bad Pixels Interpolated')
		# divider = make_axes_locatable(ax[2])
		# cax = divider.new_vertical(size="3%", pad=0.3, pack_start=True)
		# fig.add_axes(cax)
		# fig.colorbar(im3, cax=cax, orientation="horizontal", label='ADUs')
		# plt.tight_layout()
		# pdb.set_trace()
		# qp(frame_red)
		# pdb.set_trace()

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
				
	print('')
	if upload:
		print('Beginning upload process to pines.bu.edu...')
		print('Note, only PINES admins will be able to upload.')
		print('WARNING: This will overwrite data already on the PINES server!')
		print('')
		time.sleep(3)
		sftp.chdir('/data/reduced/mimir/')
		files_to_upload = np.array(natsort.natsorted(np.array([x for x in reduced_path.glob('*.fits')])))
		for i in range(len(files_to_upload)):
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
			
			print('Uploading to {}/{}, {} of {}'.format(sftp.getcwd(),nights[ind]+'/'+file.name, i+1,len(files_to_upload)))
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

	print('reduce runtime: ', round((time.time()-t1)/60,1), ' minutes.')
	print('Reduction script has completed.')
	print('')