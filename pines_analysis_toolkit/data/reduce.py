import pdb
from pathlib import Path
import os, glob
from datetime import datetime
import numpy as np
import time
from astropy.io import fits
import pdb
from pathlib import Path
from astropy.stats import sigma_clipped_stats
import matplotlib.pyplot as plt
import natsort
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.data.master_flat_chooser import master_flat_chooser
from pines_analysis_toolkit.data.master_dark_chooser import master_dark_chooser
from pines_analysis_toolkit.data.bpm_chooser import bpm_chooser
from pines_analysis_toolkit.utils.quick_plot import quick_plot as qp
from progressbar import ProgressBar
	
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
	raw_files = [Path(i) for i in natsort.natsorted(glob.glob(os.path.join(raw_path,'*.fits')))] #Natsort sorts things how you expect them to be sorted.
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
		files_to_upload = np.array(natsort.natsorted(np.array([x for x in reduced_path.glob('*.fits')])))
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