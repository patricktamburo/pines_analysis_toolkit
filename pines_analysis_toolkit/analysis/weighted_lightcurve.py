import matplotlib.pyplot as plt
colormap = plt.cm.viridis
import pdb 
from pines_analysis_toolkit.utils import pines_dir_check, short_name_creator
from pines_analysis_toolkit.analysis.night_splitter import night_splitter
from pines_analysis_toolkit.analysis.block_splitter import block_splitter
from pines_analysis_toolkit.analysis.analysis_plots import raw_flux_plot, global_raw_flux_plot, normalized_flux_plot, global_normalized_flux_plot, corr_target_plot, global_corr_target_plot, corr_all_sources_plot
from glob import glob
from natsort import natsorted
import pandas as pd
import numpy as np
import time 
import julian 
import matplotlib.dates as mdates
import os
import shutil
from copy import deepcopy
from scipy.stats import pearsonr, sigmaclip
from pathlib import Path
import os 

def weighted_lightcurve(target, phot_type='aper', convergence_threshold=1e-5, mode='night', plots=False):

    '''Authors:
		Phil Muirhead & Patrick Tamburo, Boston University, November 2020
	Purpose:
        Makes an "artificial comparison lightcurve" (ALC) using comparison star fluxes. 
	Inputs:
        target (str): the target's full 2MASS name.
        phot_type (str, optional): 'aper' for aperture photometry.
        convergence_threshold (float, optional): the threshold for maximum change for each weight between iterations before they're considered converged. By default, 1e-5 (Murray et al. 2020). 
        mode (str): either 'night' or 'global'. 'night' will make a separate lightcurve for each night of data for the target, while 'global' makes a lightcurve using all nights of data simultaneously. 
    Outputs:

	TODO:
        Make compatible with PSF photometry.
    '''
    print('\nRunning weighted_lightcurve().\n')

    if (phot_type != 'aper') and (phot_type != 'psf'):
        raise ValueError("phot_type must be either 'aper' or 'psf'!")

    #Set up paths and get photometry files. 
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    phot_path = pines_path/('Objects/'+short_name+'/'+phot_type+'_phot/')
    phot_files = natsorted(glob(str(phot_path/'*.csv')))
    analysis_path = pines_path/('Objects/'+short_name+'/analysis/')
    
    #Get source names for this observation. 
    source_detection_file = pines_path/('Objects/'+short_name+'/sources/target_and_references_source_detection.csv')
    source_df = pd.read_csv(source_detection_file)
    source_names = np.array(source_df['Name'])
    ref_names = source_names[1:]
    num_refs = len(ref_names)

    best_binned_std = 9999 #Initialize. 

    #Loop over all photometry files.
    for i in range(len(phot_files)):

        phot_file = Path(phot_files[i])
        if phot_type == 'aper':
            ap_rad = phot_file.name.split('_')[3]
        else:
            ap_rad == ''
        print('\nAperture radius: {} pixels.'.format(ap_rad))

        df = pd.read_csv(phot_file)
        df.columns = df.keys().str.strip()
        
        full_times = np.array(df['Time JD'])        

        #Split data up into individual nights. 
        if mode == 'night':
            night_inds = night_splitter(full_times)
        #Or treat all observations as a single "night".
        elif mode == 'global':
            night_inds = [list(np.arange(len(full_times)))]
        else:
            raise ValueError("mode must be either 'night' or 'global'!")

        num_nights = len(night_inds)

        #Declare some lists to hold data for all nights after processing is complete. 
        all_nights_times                 = []
        all_nights_binned_times          = []

        all_nights_raw_targ_flux         = []
        all_nights_raw_targ_err          = []
        all_nights_norm_targ_flux        = []
        all_nights_norm_targ_err         = []

        all_nights_raw_ref_flux          = []
        all_nights_raw_ref_err           = []
        all_nights_norm_ref_flux         = []
        all_nights_norm_ref_err          = []

        all_nights_alc_flux              = []
        all_nights_alc_err               = []

        all_nights_corr_targ_flux        = []
        all_nights_corr_targ_err         = []
        all_nights_binned_corr_targ_flux = []
        all_nights_binned_corr_targ_err  = []

        all_nights_corr_ref_flux         = []
        all_nights_corr_ref_err          = []
        all_nights_binned_corr_ref_flux  = []
        all_nights_binned_corr_ref_err   = []

        norm_night_weights = []

        #Loop over each night in the dataset.
        for j in range(num_nights):
            num_frames = len(night_inds[j])
            inds = np.array(night_inds[j])

            #Grab the target flux and error. 
            targ_flux = np.array(df[short_name+' Flux'][inds])
            targ_flux_err = np.array(df[short_name+' Flux Error'][inds])
            
            #Normalize the target flux and error. 
            normalization = np.sum(targ_flux/targ_flux_err**2) / np.sum(1/targ_flux_err**2)
            targ_flux_norm = targ_flux/normalization 

            #Sigmaclip the normalized target flux, in order to toss any suspect measurements. 
            vals, lower, upper = sigmaclip(targ_flux_norm, low=3, high=3)
            good_frames = np.where((targ_flux_norm > lower) & (targ_flux_norm < upper))[0]
                
            num_frames = len(good_frames)
            inds = inds[good_frames]

            #Regrab raw flux using *good* frames
            targ_flux = np.array(df[short_name+' Flux'][inds])
            targ_flux_err = np.array(df[short_name+' Flux Error'][inds])
            all_nights_raw_targ_flux.append(targ_flux)
            all_nights_raw_targ_err.append(targ_flux_err)
            
            #Recompute the normalization using only the good frames. 
            good_frame_normalization = np.sum(targ_flux/targ_flux_err**2) / np.sum(1/targ_flux_err**2)
            targ_flux_norm = targ_flux / good_frame_normalization
            targ_flux_err_norm = targ_flux_err / good_frame_normalization
            all_nights_norm_targ_flux.append(targ_flux_norm)
            all_nights_norm_targ_err.append(targ_flux_err_norm)

            #make numpy arrays for all the stuff we need
            raw_flux =       np.zeros((num_frames,num_refs)) #units of photons
            raw_err =        np.zeros((num_frames,num_refs)) #units of photons
            norm_flux =      np.zeros((num_frames,num_refs)) #raw flux normalized by a single value to be near 1.0, unitless
            norm_err =       np.zeros((num_frames,num_refs)) #norm flux err, unitless
            alc_flux =       np.zeros((num_frames,num_refs)) #special ALC to remove atm veriation, unitless, near 1.0
            alc_err =        np.zeros((num_frames,num_refs)) #err in special ALC
            corr_flux =      np.zeros((num_frames,num_refs)) #norm flux corrected to remove atm veriation, unitless, near 1.0
            corr_err =       np.zeros((num_frames,num_refs)) #corr flux err
            running_stddev = np.zeros((num_frames,num_refs)) #running std dev of corr flux

            times = np.array(full_times[inds])
            all_nights_times.append(times)

            #times = np.array([julian.from_jd(times[i], fmt='jd') for i in range(len(times))])

            #fill raw flux and err flux arrays, and normalize
            for k in np.arange(num_refs):
                r = ref_names[k]
                raw_flux[:,k] = np.array(df[r+' Flux'][inds])
                raw_err[:,k] = np.array(df[r+ ' Flux Error'][inds])

                #Normalize to weighted mean of flux values so they are near 1.
                #Weights are 1/raw_err**2.
                normalization = np.sum(raw_flux[:,k]/raw_err[:,k]**2) / np.sum(1/raw_err[:,k]**2)
                norm_flux[:,k] = raw_flux[:,k] / normalization
                norm_err[:,k] = raw_err[:,k] / normalization   
            
            all_nights_raw_ref_flux.append(raw_flux)
            all_nights_raw_ref_err.append(raw_err)
            all_nights_norm_ref_flux.append(norm_flux)
            all_nights_norm_ref_err.append(norm_err)

            #now calculate the "special" ALC for each ref star. 
            old_stddev = np.zeros(num_refs)
            new_stddev = np.zeros(num_refs)
            keep_going = True
            while keep_going:
                for k in np.arange(num_refs):
                    #indices w/o the k'th ref star
                    indices = np.where(np.arange(num_refs) != k)
                    
                    #make local arrays to store non_k norm fluxes
                    norm_flux_without_k = np.resize(norm_flux[:,indices],(num_frames,num_refs-1))
                    norm_err_without_k = np.resize(norm_err[:,indices],(num_frames,num_refs-1))

                    #ALC = weighted mean of norm flux, weights = 1/norm_err**2
                    alc_flux[:,k] = np.sum(norm_flux_without_k/norm_err_without_k**2,axis=1) / np.sum(1/norm_err_without_k**2,axis=1)
                    alc_err[:,k] = np.sqrt( 1 / np.sum(1/norm_err_without_k**2,axis=1))
                    
                    #calculate corr flux and corr err
                    corr_flux[:,k] = norm_flux[:,k] / alc_flux[:,k]
                    corr_err[:,k] = corr_flux[:,k] * np.sqrt((norm_err[:,k]/norm_flux[:,k])**2 + (alc_err[:,k]/alc_flux[:,k])**2)
                
                    #calculate stddevs
                    old_stddev[k] = new_stddev[k]
                    new_stddev[k] = np.std(corr_flux[:,k])

                    #update normalized errors
                    norm_err[:,k] = new_stddev[k]
                
                fractional_changes = np.abs(new_stddev - old_stddev)/old_stddev
                
                if np.sum(fractional_changes) < 0.000000001:
                    keep_going = False
            
            all_nights_corr_ref_flux.append(corr_flux)
            all_nights_corr_ref_err.append(corr_err)

            alc_final_flux = np.sum(norm_flux/norm_err**2,axis=1) / np.sum(1/norm_err**2,axis=1)
            alc_final_err = np.sqrt( 1 / np.sum(1/norm_err**2,axis=1))    

            norm_night_weights.append((1/norm_err[0,:])**2/np.sum((1/norm_err[0,:])**2)) #Save the NORMALIZED weight for each comparison star on this night.
            all_nights_alc_flux.append(alc_final_flux)
            all_nights_alc_err.append(alc_final_err)

            targ_flux_corr = targ_flux_norm/alc_final_flux
            targ_flux_corr_err = np.sqrt( (targ_flux_err_norm/alc_final_flux)**2 + (targ_flux_norm*alc_final_err/(alc_final_flux**2))**2 )
            all_nights_corr_targ_flux.append(targ_flux_corr)
            all_nights_corr_targ_err.append(targ_flux_corr_err)

            block_inds = block_splitter(times, [])  
            binned_times = np.zeros(len(block_inds))
            binned_flux  = np.zeros(len(block_inds))
            binned_errs  = np.zeros(len(block_inds))

            binned_ref_fluxes = np.zeros((len(block_inds), num_refs))
            binned_ref_flux_errs = np.zeros((len(block_inds), num_refs))

            for k in range(len(block_inds)):
                binned_times[k] = np.mean(times[block_inds[k]])
                binned_flux[k]  = np.mean(targ_flux_corr[block_inds[k]])
                binned_errs[k] = np.mean(targ_flux_corr_err[block_inds[k]]) / np.sqrt(len(block_inds[k])) #Uses calculated errors. 
                #binned_errs[k]  = np.std(targ_flux_corr[block_inds[k]])/np.sqrt(len(block_inds[k])) #Uses standard deviation of the block.
                
                for l in range(num_refs):
                    binned_ref_fluxes[k,l] = np.mean(corr_flux[:,l][block_inds[k]])
                    binned_ref_flux_errs[k,l] = np.mean(corr_err[:,l][block_inds[k]]) / np.sqrt(len(block_inds[k])) #Uses calculated errors.
                    #binned_ref_flux_errs[k,l] = np.std(corr_flux[:,l][block_inds[k]])/np.sqrt(len(block_inds[k])) #Uses standard deviation of the block.
            
            all_nights_binned_times.append(binned_times)
            all_nights_binned_corr_targ_flux.append(binned_flux)
            all_nights_binned_corr_targ_err.append(binned_errs)
            all_nights_binned_corr_ref_flux.append(binned_ref_fluxes)
            all_nights_binned_corr_ref_err.append(binned_ref_flux_errs)

        #Unpack the binned target flux into a list for ease of calculating stddev.    
        all_bin_flux = []
        for kk in range(len(all_nights_binned_corr_targ_flux)):
            all_bin_flux.extend(all_nights_binned_corr_targ_flux[kk])

        binned_std = np.std(all_bin_flux)
        print('Stddev of binned corrected target flux: {:1.4f}'.format(binned_std))
        if binned_std < best_binned_std:
            print('New best lightcurve found!')
            best_binned_std = binned_std
            best_ap = ap_rad

        #Output the weighted lc data to a csv. 
        time_save = times
        flux_save = targ_flux_corr 
        flux_err_save = targ_flux_corr_err
        output_dict = {'Time':time_save, 'Flux':flux_save, 'Flux Error':flux_err_save}
        output_df = pd.DataFrame(data=output_dict)
        if phot_type == 'aper':
            if not os.path.exists(analysis_path/('aper_phot_analysis/'+ap_rad)):
                os.mkdir(analysis_path/('aper_phot_analysis/'+ap_rad))
            output_filename = analysis_path/('aper_phot_analysis/'+ap_rad+'/'+short_name+'_weighted_lc_aper_phot_'+ap_rad+'_pix.csv')
            print('Saving to {}.'.format(output_filename))
            output_df.to_csv(output_filename)
        elif phot_type == 'psf':
            print("ERROR: Need to create flux output for PSF photometry.")     

        #TODO: Would be nice to wrap all this up in a dictionary that gets passed to the plotting programs. 
        #Plot the raw fluxes, normalized fluxes, corrected target flux, and corrected reference star fluxes. 
        if plots:
            if mode == 'night':
                raw_flux_plot(all_nights_times, all_nights_raw_targ_flux, all_nights_raw_targ_err, all_nights_raw_ref_flux, all_nights_raw_ref_err, short_name, analysis_path, phot_type, ap_rad)   
                normalized_flux_plot(all_nights_times, all_nights_norm_targ_flux, all_nights_norm_targ_err, all_nights_norm_ref_flux, all_nights_norm_ref_err, all_nights_alc_flux, all_nights_alc_err, short_name, analysis_path, phot_type, ap_rad)
                corr_target_plot(all_nights_times, all_nights_corr_targ_flux, all_nights_binned_times, all_nights_binned_corr_targ_flux, all_nights_binned_corr_targ_err, short_name, analysis_path, phot_type, ap_rad)
                corr_all_sources_plot(all_nights_times, all_nights_corr_targ_flux, all_nights_binned_times, all_nights_binned_corr_targ_flux, all_nights_binned_corr_targ_err, all_nights_corr_ref_flux, all_nights_binned_corr_ref_flux, all_nights_binned_corr_ref_err, short_name, analysis_path, phot_type, ap_rad, num_refs, num_nights, norm_night_weights)
            elif mode == 'global':
                global_raw_flux_plot(all_nights_times, all_nights_raw_targ_flux, all_nights_raw_targ_err, all_nights_raw_ref_flux, all_nights_raw_ref_err, short_name, analysis_path, phot_type, ap_rad)
                global_normalized_flux_plot(all_nights_times, all_nights_norm_targ_flux, all_nights_norm_targ_err, all_nights_norm_ref_flux, all_nights_norm_ref_err, all_nights_alc_flux, all_nights_alc_err, short_name, analysis_path, phot_type, ap_rad)
                global_corr_target_plot(all_nights_times, all_nights_corr_targ_flux, all_nights_binned_times, all_nights_binned_corr_targ_flux, all_nights_binned_corr_targ_err, short_name, analysis_path, phot_type, ap_rad)
                
    #Write the best aperture to a text file.
    print('\nBest aperture: {}.'.format(best_ap))
    best_ap_output_path = analysis_path/('optimal_aperture.txt')
    with open(best_ap_output_path, 'w') as f:
        f.write(best_ap)
    
    return