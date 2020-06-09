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


'''Authors:
		Paul Dalba, Boston University, February 2017
		Patrick Tamburo, Boston University, July 2019 and June 2020
	Purpose:
		Reduces raw PINES science images for a specified target and saves the reduced images. 
	Inputs:
		target_name (str): the target's full 2MASS name, i.e. '2MASS J01234567+0123456' 
		flat_type (str, optional): either 'dome' or 'sky'. By default, set to dome. We haven't seen better performance with sky flats. 
	Outputs:
		None
	TODO:
		None
'''
	

def reduce(target_name, flat_type='dome'):

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
	raw_path = pines_path+'Objects/'+short_name+'/raw/'
	raw_files = natsort.natsorted(glob.glob(raw_path+'*.fits')) #Natsort sorts things how you expect them to be sorted.
	dark_path = pines_path+'Calibrations/Darks/'
	flats_path = pines_path+'Calibrations/Flats/'
	if flat_type == 'sky':
		reduced_path = pines_path+'Objects/'+short_name+'/skyflat_reduced/'
	elif flat_type == 'dome':
		reduced_path = pines_path+'Objects/'+short_name+'/domeflat_reduced/'
		
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
		pdb.set_trace()
		master_dark = fits.open(dark_path+'master_dark_'+str(exptime)+'.fits')[0].data 
		if flat_type == 'dome': #Read in the appropriate flat. 
			possible_flats = glob.glob(flats_path+'Domeflats/'+band+'/Master flats/*.fits')
			flat_dates = [datetime.strptime(possible_flats[i].split('_')[-1].split('.')[0],'%Y%m%d') for i in range(len(possible_flats))]
			file_date = datetime.strptime(header['DATE-OBS'].split('T')[0].replace('-',''),'%Y%m%d')
			closest_flat =  min(flat_dates, key=lambda d: abs(d-file_date)).strftime('%Y%m%d')
			master_flat = fits.open(flats_path+'Domeflats/'+band+'/Master flats/master_flat_'+band+'_'+closest_flat+'.fits')[0].data
			print('Using master_flat_'+band+'_'+closest_flat+'.fits')
		elif flat_type == 'sky':
			print('ERROR: Still have to add skyflat functionality.')
			pdb.set_trace()
		
		frame_corr1 = frame_raw - master_dark 
		
		#Flat field correction
		frame_red = frame_corr1/master_flat

		#Save the new flat file- remove the previous file if it does exist
		fits.writeto(target_filename, frame_red, header)

		print('')
		print('Reducing '+raw_files[i].split('/')[-1]+': '+str(i+1)+' of '+str(size(raw_files)),', band = ',band,' Exptime = ',exptime)
		
	#End of script indicator
	print('')
	print('reduce runtime: ', np.round((time.time()-t1)/60,1), ' minutes.')
	print('Reduction script has completed.')
	print('')