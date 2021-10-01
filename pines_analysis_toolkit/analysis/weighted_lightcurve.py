import matplotlib.pyplot as plt
colormap = plt.cm.viridis
import pdb 
from pines_analysis_toolkit.utils import pines_dir_check, short_name_creator
from pines_analysis_toolkit.analysis.night_splitter import night_splitter
from pines_analysis_toolkit.analysis.block_splitter import block_splitter
from pines_analysis_toolkit.analysis.block_binner import block_binner
from pines_analysis_toolkit.analysis.analysis_plots import *
from pines_analysis_toolkit.analysis.regression import *
from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader
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
import fileinput
import time   

def weighted_lightcurve(short_name, phot_type='aper', convergence_threshold=1e-9, mode='night', n_sig_refs=5, sigma_clip_threshold=4., max_iterations=1000, use_pwv=False, red_stars_only=False, blue_stars_only=False, high_correlation_refs=False, force_output_path=''):
    """Creates a corrected light curve for the target using a weighted mean of reference star fluxes. 

    :param short_name: short name of the target
    :type short_name: str
    :param phot_type: photometry type, 'aper' or 'psf', defaults to 'aper'
    :type phot_type: str, optional
    :param convergence_threshold: threshold in reference star weighting loop to consider weights converged, defaults to 1e-9
    :type convergence_threshold: float, optional
    :param mode: 'night' or 'global' normalization mode, defaults to 'night'
    :type mode: str, optional
    :param n_sig_refs: n brightest reference stars whose average standard deviation is evaluated to choose the best aperture, defaults to 5
    :type n_sig_refs: int, optional
    :param sigma_clip_threshold: threshold for sigma clipping of target light curve, defaults to 4.0
    :type sigma_clip_threshold: float, optional
    :param max_iterations: maximum number of iterations allowed in reference star weighting loop, defaults to 1000
    :type max_iterations: int, optional
    :param use_pwv: whether or not to use precipitable water vapor data in the regression, defaults to False
    :type use_pwv: bool, optional
    :param red_stars_only: whether or not to use only the reddest stars in the field as references, defaults to False
    :type red_stars_only: bool, optional
    :param blue_stars_only: whether or not to use only the bluest stars in the field as references, defaults to False
    :type blue_stars_only: bool, optional
    :param high_correlation_refs: whether or not to use only the references that are most highly correlated with the raw target flux, defaults to False
    :type high_correlation_refs: bool, optional
    :param force_output_path: user-chosen path to use if you do not want to use the default ~/Documents/PINES_analysis_toolkit directory for analysis, defaults to ''
    :type force_output_path: str, optional
    """
    print('\nRunning weighted_lightcurve() in {} normalization mode.'.format(mode))

    if (phot_type != 'aper') and (phot_type != 'psf'):
        raise ValueError("phot_type must be either 'aper' or 'psf'!")
    
    def optimal_aperture_output():
        '''
            PURPOSE: writes the aperture type that minimizes the eference star scatter metric to a file 'optimal_aperture.txt' in the object's analysis directory.
        '''
        best_ap_output_path = analysis_path/('optimal_aperture.txt')
        if not os.path.exists(best_ap_output_path):
            with open(best_ap_output_path, 'w') as f:
                f.write('Nightly-normalized:  \n')
                f.write('Globally-normalized: ')

        #Write the best aperture to a text file.
        if mode == 'night':
            print('\nBest aperture for nightly-normalized photometry: {}.'.format(best_ap))
            with open(best_ap_output_path, 'r') as f:
                lines = f.readlines()
            
            lines[0] = 'Nightly-normalized:  '+best_ap+'\n'
            with open(best_ap_output_path, 'w') as f:
                for line in lines:
                    f.write(line)

        elif mode == 'global':
            print('\nBest aperture for globally-normalized photometry: {}.'.format(best_ap))
            with open(best_ap_output_path, 'r') as f:
                lines = f.readlines()
                    
            lines[1] = 'Globally-normalized: '+best_ap
            with open(best_ap_output_path, 'w') as f:
                for line in lines:
                    f.write(line)

    def optimizing_metric_caclulator(all_nights_binned_corr_ref_flux, n_sig_refs):
        all_ref_bin_avg_std = [[] for x in range(len(all_nights_binned_corr_ref_flux[0][0]))]
        for kk in range(num_nights):
            for jj in range(num_refs):
                all_ref_bin_avg_std[jj].extend([np.nanmean(np.nanstd(all_nights_binned_corr_ref_flux[kk][:,jj]))])

        #Use the average standard deviation of the binned data for the brightest n_sig_refs reference stars to determine the optimal aperture radius, fixed or variable. 
        #We use multiple reference stars (default = 5) to calculate this statistic, in case of any variability in the reference stars. 
        #We want to *minimize* binned_std to determine the best radius aperture. 
        return np.mean(all_ref_bin_avg_std[0:n_sig_refs])

    def output_filename_generator(analysis_path, ap_rad, mode):
        if phot_type == 'aper':
            if not os.path.exists(analysis_path/('aper_phot_analysis/'+ap_rad)):
                os.mkdir(analysis_path/('aper_phot_analysis/'+ap_rad))
            if 'fixed' in ap_rad:
                rad = ap_rad.split('_')[0]
                if mode == 'night':
                    output_filename = analysis_path/('aper_phot_analysis/'+ap_rad+'/'+short_name.replace(' ','')+'_fixed_aper_phot_r='+rad+'_nightly_weighted_lc.csv')
                elif mode == 'global':
                    output_filename = analysis_path/('aper_phot_analysis/'+ap_rad+'/'+short_name.replace(' ','')+'_fixed_aper_phot_r='+rad+'_global_weighted_lc.csv')

            elif 'variable' in ap_rad:
                factor = ap_rad.split('_')[0]
                if mode == 'night':
                    output_filename = analysis_path/('aper_phot_analysis/'+ap_rad+'/'+short_name.replace(' ','')+'_variable_aper_phot_f='+factor+'_nightly_weighted_lc.csv')
                elif mode == 'global':
                    output_filename = analysis_path/('aper_phot_analysis/'+ap_rad+'/'+short_name.replace(' ','')+'_variable_aper_phot_f='+factor+'_global_weighted_lc.csv')
        else:
            raise RuntimeError('Need to implement PSF photometry!')
        return output_filename

    np.seterr(divide='ignore') #Ignore divide-by-zero errors. 

    #Set up paths and get photometry files. 
    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()

    phot_path = pines_path/('Objects/'+short_name+'/'+phot_type+'_phot/')
    phot_files = natsorted(glob(str(phot_path/'*.csv')))
    analysis_path = pines_path/('Objects/'+short_name+'/analysis/')
    
    #Get source names for this observation. 
    source_detection_file = pines_path/('Objects/'+short_name+'/sources/target_and_references_source_detection.csv')
    source_df = pd.read_csv(source_detection_file)

    if red_stars_only and blue_stars_only:
        raise ValueError('red_stars_only and blue_stars_only cannot both be True!')

    #If red_stars_only, limit the source_df to the red stars.
    if red_stars_only:
        print('Using red reference stars only!')
        target_entry = source_df.loc[[0]]

        #Grab only the reddest sources from the sources DataFrame.
        source_df = source_df[(source_df['M_G'] < 14) & (source_df['M_G'] > 7.5)]

        #Create the new source_df with the target in the first row. 
        source_df = pd.concat([target_entry, source_df])

        #Reset the indices of the DataFrame.s
        source_df.reset_index(inplace=True, drop=True)

    if blue_stars_only:
        print('Using blue reference stars only!')
        target_entry = source_df.loc[[0]]

        #Grab only the reddest sources from the sources DataFrame.
        source_df = source_df[(source_df['M_G'] < 7.5)]

        #Create the new source_df with the target in the first row. 
        source_df = pd.concat([target_entry, source_df])

        #Reset the indices of the DataFrame.s
        source_df.reset_index(inplace=True, drop=True)
    
    source_names = np.array(source_df['Name'])
    ref_names = source_names[1:]
    num_refs = len(ref_names)

    #Get centroids for regression. 
    centroid_path = pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv')
    centroid_df = pines_log_reader(centroid_path)
    full_centroid_x = np.array(centroid_df[source_names[0]+' Image X'], dtype='float')
    full_centroid_y = np.array(centroid_df[source_names[0]+' Image Y'], dtype='float')

    if use_pwv:
        #Get PWV for regression. 
        pwv_path = pines_path/('Objects/'+short_name+'/pwv/'+short_name+'_fyodor_pwv.csv')
        pwv_df = pines_log_reader(pwv_path)
        full_pwv = np.array(pwv_df['PWV'])

    best_metric = 9999 #Initialize the metric that we'll use to choose the best lightcurve.

    #Loop over all photometry files.
    for i in range(len(phot_files)):
        phot_file = Path(phot_files[i])
        if phot_type == 'aper':
            if 'variable' in phot_file.name:
                ap_rad = phot_file.name.split('aper_phot_')[1].split('_')[0]+'_variable'
                print('\nRestoring variable aperture photometry output, seeing multiplicative factor = {}.'.format(phot_file.name.split('_')[4]))

            elif 'fixed' in phot_file.name:
                ap_rad = phot_file.name.split('aper_phot_')[1].split('_')[0]+'_fixed'
                print('\nRestoring fixed aperture photometry output, pixel radius = {}.'.format(ap_rad.split('_')[0]))
        else:
            ap_rad == ''

        df = pd.read_csv(phot_file)
        df.columns = df.keys().str.strip()
        
        full_times = np.array(df['Time BJD TDB'])        
        full_file_list = np.array(df['Filename'])
        full_airmass = np.array(df['Airmass']) #Get airmass for regression.
        full_background = np.array(df[source_names[0]+' Background'], dtype='float')
        sigma_clip_flags = np.zeros(len(full_times), dtype='int') #Track any frames identified as bad in sigma clipping. 

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
        all_nights_file_list             = []
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
            num_frames = len(inds)
            airmass = full_airmass[inds]
            centroid_x = full_centroid_x[inds]
            centroid_y = full_centroid_y[inds]
            background = full_background[inds]
            if use_pwv:
                pwv = full_pwv[inds]     
                        
            times = np.array(full_times[inds])
            all_nights_times.append(times)
            date = julian.from_jd(times[0])
            date_str = 'UT '+date.strftime('%b %d, %Y')

            #Grab the target flux and error. 
            targ_flux = np.array(df[source_names[0]+' Flux'][inds], dtype='float')
            targ_flux_err = np.array(df[source_names[0]+' Flux Error'][inds], dtype='float')
            
            #Normalize the target flux and error. 
            normalization = np.nansum(targ_flux/targ_flux_err**2) / np.nansum(1/targ_flux_err**2)
            targ_flux_norm = targ_flux/normalization 

            #Sigmaclip the normalized target flux, in order to toss any suspect measurements. 
            vals, lower, upper = sigmaclip(targ_flux_norm, low=sigma_clip_threshold, high=sigma_clip_threshold)
            good_frames = np.where((targ_flux_norm > lower) & (targ_flux_norm < upper))[0]
            bad_frames = np.where((targ_flux_norm < lower) | (targ_flux_norm > upper))[0]
                
            #inds = inds[good_frames]
            bad_inds = inds[bad_frames]
            sigma_clip_flags[bad_inds] = 1 #Flag the frame as bad (1) for the output. 

            #Regrab raw flux using *good* frames. 
            targ_flux[bad_frames] = np.nan
            targ_flux_err[bad_frames] = np.nan

            #targ_flux = np.array(df[short_name+' Flux'][inds])
            #targ_flux_err = np.array(df[short_name+' Flux Error'][inds])
            all_nights_raw_targ_flux.append(targ_flux)
            all_nights_raw_targ_err.append(targ_flux_err)
          
            #Recompute the normalization using only the good frames. 
            ##good_frame_normalization = np.nansum(targ_flux/targ_flux_err**2) / np.nansum(1/targ_flux_err**2)
            good_frame_normalization = np.nanmean(targ_flux)
            targ_flux_norm = targ_flux / good_frame_normalization
            targ_flux_err_norm = targ_flux_err / good_frame_normalization
            all_nights_norm_targ_flux.append(targ_flux_norm)
            all_nights_norm_targ_err.append(targ_flux_err_norm)

            #make numpy arrays for all the stuff we need
            raw_flux =       np.zeros((num_frames,num_refs)) #units of photons
            raw_err =        np.zeros((num_frames,num_refs)) #units of photons
            norm_flux =      np.zeros((num_frames,num_refs)) #raw flux normalized by a single value to be near 1.0, unitless
            norm_err =       np.zeros((num_frames,num_refs)) #norm flux err, unitless
            

            try:    
                block_inds = block_splitter(times) #Get indices of each block for binning/plotting purposes.
            except:
                pdb.set_trace()
                
            files = np.array(full_file_list[inds])
            all_nights_file_list.append(files)

            #Write out bad frames, if any. 
            if len(bad_frames) > 0:
                print('Sigma clipping identified {} bad frames:'.format(len(bad_frames)))
                for f in range(len(bad_frames)):
                    print('     '+full_file_list[bad_inds[f]])

            #fill raw flux and err flux arrays, and normalize
            for k in np.arange(num_refs):
                r = ref_names[k]
                raw_flux[:,k] = np.array(df[r+' Flux'][inds], dtype='float')
                raw_err[:,k] = np.array(df[r+ ' Flux Error'][inds], dtype='float')
                raw_flux[bad_frames, k] = np.nan #Overwrite any sigma-clipped frames with NaNs. 
                raw_err[bad_frames, k] = np.nan
                
                #Normalize to weighted mean of flux values so they are near 1.
                #Weights are 1/raw_err**2.
                normalization = np.nansum(raw_flux[:,k]/raw_err[:,k]**2) / np.nansum(1/raw_err[:,k]**2) #ORIGINAL VERSION!

                norm_flux[:,k] = raw_flux[:,k] / normalization
                norm_err[:,k] = raw_err[:,k] / normalization   

            all_nights_raw_ref_flux.append(raw_flux)
            all_nights_raw_ref_err.append(raw_err)
            all_nights_norm_ref_flux.append(norm_flux)
            all_nights_norm_ref_err.append(norm_err)

            if high_correlation_refs:
                selection_type = 'worst'
                n_to_take = 4

                #Measure the Pearson correlation coefficient of each reference normalized flux with the target's normalized flux.
                correlations = np.zeros(num_refs)
                for k in range(num_refs):
                    non_nan_inds = np.where(~np.isnan(norm_flux[:,k]) | ~np.isnan(targ_flux_norm))[0]
                    corr, sig = pearsonr(targ_flux_norm[non_nan_inds], norm_flux[:,k][non_nan_inds])
                    correlations[k] = corr

                #Sort the references by those with the highest correlation with the target, and take the best handful of them. 
                sorted_corr_inds = np.argsort(correlations)[::-1]

                #Can choose the n_to_take most correlated references...
                if selection_type == 'best':
                    best_ref_inds = sorted_corr_inds[0:n_to_take]
                #Or the least correlated, for testing purposes...
                elif selection_type == 'worst':
                    best_ref_inds = sorted_corr_inds[n_to_take:]

                #Change the norm_flux array to only contain the flux from the most highly correlated references.
                norm_flux = norm_flux[:,best_ref_inds]
                norm_err  = norm_err[:,best_ref_inds]
                num_refs = len(norm_flux[0,:])
                ref_names = ref_names[best_ref_inds]

            #Define some more arrays. 
            alc_flux =       np.zeros((num_frames,num_refs)) #special ALC to remove atm veriation, unitless, near 1.0
            alc_err =        np.zeros((num_frames,num_refs)) #err in special ALC
            corr_flux =      np.zeros((num_frames,num_refs)) #norm flux corrected to remove atm veriation, unitless, near 1.0
            corr_err =       np.zeros((num_frames,num_refs)) #corr flux err
            regressed_corr_flux = np.zeros((num_frames,num_refs)) #corrected flux that has been run through the linear regression procedure

            #Now calculate the "special" ALC for each ref star, use that to correct the ref star's flux, and run the corrected flux through the linear regression model.
            #Regressed corrected flux is used to measure the standard deviation of the lightcurve, which is used to weight it in the target ALC. 
            old_stddev = np.zeros(num_refs)
            new_stddev = np.zeros(num_refs)
            keep_going = True
            iteration = 0 
            while keep_going:
                iteration += 1

                for k in np.arange(num_refs):
                    #indices w/o the k'th ref star
                    indices = np.where(np.arange(num_refs) != k)
                    
                    #make local arrays to store non_k norm fluxes
                    norm_flux_without_k = np.resize(norm_flux[:,indices],(num_frames,num_refs-1))
                    norm_err_without_k = np.resize(norm_err[:,indices],(num_frames,num_refs-1))

                    #ALC = weighted mean of norm flux, weights = 1/norm_err**2
                    alc_flux[:,k] = np.sum(norm_flux_without_k/norm_err_without_k**2,axis=1) / np.sum(1/norm_err_without_k**2,axis=1)
                    alc_err[:,k] = np.sqrt( 1 / np.sum(1/norm_err_without_k**2,axis=1))
                          
                    #Renormalize the ALC and normalized reference flux so that their mean is EXACTLY 1. 
                    alc_renorm = np.nanmean(alc_flux[:,k])
                    alc_flux[:,k] = alc_flux[:,k] / alc_renorm
                    norm_renorm = np.nanmean(norm_flux[:,k])
                    norm_flux[:,k] = norm_flux[:,k] / norm_renorm

                    #calculate corr flux and corr err
                    corr_flux[:,k] = norm_flux[:,k] / alc_flux[:,k]
                    corr_err[:,k] = corr_flux[:,k] * np.sqrt((norm_err[:,k]/norm_flux[:,k])**2 + (alc_err[:,k]/alc_flux[:,k])**2)

                    #Perform a regression of the corrected flux against various parameters. 
                    if use_pwv:
                        regression_dict = {'airmass':airmass, 'centroid_x':centroid_x, 'centroid_y':centroid_y, 'pwv':pwv}
                    else:
                        regression_dict = {'airmass':airmass, 'centroid_x':centroid_x, 'centroid_y':centroid_y}

                    regressed_corr_flux[:,k] = regression(corr_flux[:,k], regression_dict)

                    #Calculate stddevs from the regressed corrected flux.
                    old_stddev[k] = new_stddev[k]
                    new_stddev[k] = np.nanstd(regressed_corr_flux[:,k])

                    #update normalized errors
                    norm_err[:,k] = new_stddev[k]

                fractional_changes = np.abs(new_stddev - old_stddev)/old_stddev
                
                #Once fractional changes stabilize, break out of the loop. 
                if np.sum(fractional_changes) < convergence_threshold:
                    keep_going = False

                if iteration > max_iterations:
                    print('Iterations in reference star weighting loop exceeded max_iterations.')
                    break
                
            
            all_nights_corr_ref_flux.append(regressed_corr_flux)
            all_nights_corr_ref_err.append(corr_err)

            alc_final_flux = np.sum(norm_flux/norm_err**2,axis=1) / np.sum(1/norm_err**2,axis=1)
            alc_final_err = np.sqrt( 1 / np.sum(1/norm_err**2,axis=1))    

            norm_night_weights.append((1/norm_err[0,:])**2/np.sum((1/norm_err[0,:])**2)) #Save the NORMALIZED weight for each comparison star on this night.
            all_nights_alc_flux.append(alc_final_flux)
            all_nights_alc_err.append(alc_final_err)

            #Correct the target flux with the master ALC and renormalize to 1. 
            targ_flux_corr = targ_flux_norm/alc_final_flux
            targ_flux_corr_err = np.sqrt( (targ_flux_err_norm/alc_final_flux)**2 + (targ_flux_norm*alc_final_err/(alc_final_flux**2))**2)
            
            if use_pwv:
                regression_dict = {'airmass':airmass, 'centroid_x':centroid_x, 'centroid_y':centroid_y,'pwv':pwv}
            else: 
                regression_dict = {'airmass':airmass, 'centroid_x':centroid_x, 'centroid_y':centroid_y}
            
            regressed_targ_flux_corr = regression(targ_flux_corr, regression_dict, verbose=True)

            
            #regressed_targ_flux_corr = leave_one_out_regression(times, targ_flux_corr, targ_flux_corr_err, regression_dict, verbose=False)
            #pdb.set_trace()

            all_nights_corr_targ_flux.append(regressed_targ_flux_corr)
            all_nights_corr_targ_err.append(targ_flux_corr_err)

            binned_times = np.zeros(len(block_inds))
            binned_flux  = np.zeros(len(block_inds))
            binned_errs  = np.zeros(len(block_inds))

            binned_ref_fluxes = np.zeros((len(block_inds), num_refs))
            binned_ref_flux_errs = np.zeros((len(block_inds), num_refs))

            for k in range(len(block_inds)):
                binned_times[k] = np.nanmean(times[block_inds[k]])
                binned_flux[k]  = np.nanmean(regressed_targ_flux_corr[block_inds[k]])
                #binned_errs[k] = np.nanmean(targ_flux_corr_err[block_inds[k]]) / np.sqrt(len(block_inds[k])) #Uses calculated errors. 
                binned_errs[k]  = np.nanstd(targ_flux_corr[block_inds[k]])/np.sqrt(len(block_inds[k])) #Uses standard deviation of the block.
                
                for l in range(num_refs):
                    binned_ref_fluxes[k,l] = np.nanmean(regressed_corr_flux[:,l][block_inds[k]])
                   #binned_ref_flux_errs[k,l] = np.nanmean(corr_err[:,l][block_inds[k]]) / np.sqrt(len(block_inds[k])) #Uses calculated errors.
                    binned_ref_flux_errs[k,l] = np.std(corr_flux[:,l][block_inds[k]])/np.sqrt(len(block_inds[k])) #Uses standard deviation of the block.

            all_nights_binned_times.append(binned_times)
            all_nights_binned_corr_targ_flux.append(binned_flux)
            all_nights_binned_corr_targ_err.append(binned_errs)
            all_nights_binned_corr_ref_flux.append(binned_ref_fluxes)
            all_nights_binned_corr_ref_err.append(binned_ref_flux_errs)

        #Calculate our metric for choosing our best lightcurve. 
        optimizing_metric = optimizing_metric_caclulator(all_nights_binned_corr_ref_flux, n_sig_refs)

        print('Average stddev of binned corrected reference flux for brightest {} references: {:1.5f}'.format(n_sig_refs, optimizing_metric))
        if optimizing_metric < best_metric:
            print('New best lightcurve found!')
            best_metric = optimizing_metric
            best_ap = ap_rad

        #Write out the best aperture 
        optimal_aperture_output()
        
        #Output the weighted lc data to a csv. 
        #Have to unpack things that were saved into night arrays. 
        
        file_list_save = full_file_list
        time_save = full_times
        sigma_clip_flag_save = sigma_clip_flags

        targ_flux_corr_save = []
        targ_flux_corr_err_save = []
        targ_flux_norm_save = []
        targ_err_norm_save = []
        alc_flux_save = []
        alc_err_save = []
        ref_flux_corr_save = [[] for x in range(num_refs)]
        ref_flux_corr_err_save = [[] for x in range(num_refs)]
        night_weight_arrays = [[] for x in range(num_refs)]
        for j in range(num_nights):
            targ_flux_corr_save.extend(all_nights_corr_targ_flux[j])
            targ_flux_corr_err_save.extend(all_nights_corr_targ_err[j])
            targ_flux_norm_save.extend(all_nights_norm_targ_flux[j])
            targ_err_norm_save.extend(all_nights_norm_targ_err[j])
            alc_flux_save.extend(all_nights_alc_flux[j])
            alc_err_save.extend(all_nights_alc_err[j])
            for k in range(num_refs):
                ref_flux_corr_save[k].extend(all_nights_corr_ref_flux[j][:,k])
                ref_flux_corr_err_save[k].extend(all_nights_corr_ref_err[j][:,k])
                night_weight_arrays[k].extend(np.zeros(len(night_inds[j]))+norm_night_weights[j][k])
        
        #Create the output dictionary
        output_dict = {'Filename':file_list_save, 'Time BJD TDB':time_save, 'Sigma Clip Flag':sigma_clip_flag_save, short_name+' Normalized Flux':targ_flux_norm_save, short_name+' Normalized Flux Error':targ_err_norm_save, short_name+' Corrected Flux':targ_flux_corr_save, short_name+' Corrected Flux Error':targ_flux_corr_err_save, 'ALC Flux':alc_flux_save, 'ALC Flux Error':alc_err_save}
        for j in range(num_refs):
            ref = ref_names[j]
            output_dict[ref+' Corrected Flux'] = ref_flux_corr_save[j]
            output_dict[ref+' Corrected Flux Error'] = ref_flux_corr_err_save[j]
        
        #Convert to a dataframe
        output_df = pd.DataFrame(data=output_dict)
        output_filename = output_filename_generator(analysis_path, ap_rad, mode)
        print('Saving weighted_lightcurve output to {}.'.format(output_filename))

        #Write out to a csv file.
        with open(output_filename, 'w') as f:
            for j in range(len(output_df)):
                #Write the header line. 
                if j == 0:
                    f.write('{:<21s}, '.format('Filename'))
                    f.write('{:<15s}, '.format('Time BJD TDB'))
                    f.write('{:<15s}, '.format('Sigma Clip Flag'))
                    f.write('{:<31s}, {:<31s}, {:<30s}, {:<30s}, {:<9s}, {:<15s}, '.format(short_name+' Normalized Flux',short_name+' Normalized Flux Error',short_name+' Corrected Flux', short_name+' Corrected Flux Error', 'ALC Flux','ALC Flux Error'))

                    for i in range(len(ref_names)):
                        n = ref_names[i]
                        if i != num_refs - 1:
                            f.write('{:<22s}, {:<30s}, {:<36s}, '.format(n+' ALC Weight', n+' Corrected Flux', n+' Corrected Flux Error'))
                        else:
                            f.write('{:<22s}, {:<30s}, {:<36s}\n'.format(n+' ALC Weight', n+' Corrected Flux', n+' Corrected Flux Error'))
                
                #Write in the data lines.
                f.write('{:<21s}, '.format(output_df['Filename'][j]))
                f.write('{:<15.7f}, '.format(output_df['Time BJD TDB'][j]))
                f.write('{:<15d}, '.format(output_df['Sigma Clip Flag'][j]))
                f.write('{:<31.6f}, {:<31.6f}, {:<30.6f}, {:<36.6f}, {:<9.6f}, {:<15.6f}, '.format(output_df[short_name+' Normalized Flux'][j],output_df[short_name+' Normalized Flux Error'][j],output_df[short_name+' Corrected Flux'][j], output_df[short_name+' Corrected Flux Error'][j], output_df['ALC Flux'][j], output_df['ALC Flux Error'][j]))
                
                for i in range(num_refs):
                    n = ref_names[i]                  
                    if i != num_refs - 1:
                        format_string = '{:<22.3f}, {:<30.6f}, {:<36.6f}, '
                    else:
                        format_string = '{:<22.3f}, {:<30.6f}, {:<36.6f}\n'
                    f.write(format_string.format(night_weight_arrays[i][j],output_df[n+' Corrected Flux'][j], output_df[n+' Corrected Flux Error'][j]))
    print('\nBest aperture: {}.\n'.format(best_ap))
    time.sleep(2)
    return