import matplotlib.pyplot as plt
colormap = plt.cm.viridis
import pdb 
from pines_analysis_toolkit.utils import pines_dir_check, short_name_creator
from pines_analysis_toolkit.analysis.night_splitter import night_splitter
from pines_analysis_toolkit.analysis.block_splitter import block_splitter
from pines_analysis_toolkit.analysis.pines_lc_plot import pines_lc_plot
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
from scipy.stats import pearsonr

def pines_alc(target, phot_type='aper', convergence_threshold=1e-5, mode='night'):

    '''Authors:
		Phil Muirhead & Patrick Tamburo, Boston University, November 2020
	Purpose:
        Makes an "artificial comparison lightcurve" (ALC) using comparison star fluxes. 
	Inputs:
        target (str): the target's full 2MASS name.
        phot_type (str, optional): 'aper' for aperture photometry.
        convergence_threshold (float, optional): the threshold for maximum change for each weight between iterations before they're considered converged. By default, 1e-5 (Murray et al. 2020). 
        mode (str): either 'night' or 'global. 'night' will make a separate lightcurve for each night of data for the target, while 'global' makes a lightcurve using all nights of data simultaneously. 
    Outputs:

	TODO:
        Make compatible with PSF photometry.
        Add 'global' performance. 
    '''

    #Set up paths and get photometry files. 
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    phot_path = pines_path/('Objects/'+short_name+'/'+phot_type+'_phot/')
    phot_files = natsorted(glob(str(phot_path/'*.csv')))
    analysis_path = pines_path/('Objects/'+short_name+'/analysis/')

    #Loop over all photometry files.
    for i in range(len(phot_files)):

        phot_file = phot_files[i]
        df = pd.read_csv(phot_file)
        df.columns = df.keys().str.strip()
        
        full_times = np.array(df['Time JD'])

        source_names = natsorted(list(set([i.split(' ')[0]+' '+i.split(' ')[1] for i in df.keys() if (i[0] == '2') or (i[0] == 'R')])))
        ref_names = source_names[1:]
        num_refs = len(ref_names)

        #Split data up into individual nights. 
        night_inds = night_splitter(full_times)
        num_nights = len(night_inds)

        #Loop over each night in the dataset.
        for j in range(num_nights):
            num_frames = len(night_inds[j])
            inds = night_inds[j]
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
            ref_x_pos =      np.zeros((num_frames,num_refs)) #The x pixel position of reference stars on the chip
            ref_y_pos =      np.zeros((num_frames,num_refs)) #The y pixel position of reference stars on the chip

            times = np.array(full_times[inds])
            #times = np.array([julian.from_jd(times[i], fmt='jd') for i in range(len(times))])

            targ_flux = np.array(df[short_name+' Flux'][inds])
            targ_flux_err = np.array(df[short_name+' Flux Error'][inds])
            normalization = np.sum(targ_flux/targ_flux_err**2) / np.sum(1/targ_flux_err**2)
            targ_flux_norm = targ_flux/normalization 
            targ_flux_err_norm = targ_flux_err / normalization

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
            
            # #Plot normalized reference star fluxes. 
            # fig, ax = plt.subplots(1,1,figsize=(14,5))
            # ax.plot(times, norm_flux, '.')
            # ax.errorbar(times, targ_flux_norm, targ_flux_err_norm, marker='o', linestyle='')
            # ax.set_title('Normalized Reference Star Fluxes', fontsize=16)
            # ax.tick_params(labelsize=12)
            # ax.set_ylabel('Normalized Flux', fontsize=14)
            # ax.set_xlabel('Time (JD$_{UTC}$)', fontsize=14)
            # plt.tight_layout()
            # plt.savefig(analysis_path/'norm_refs.png')


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
                
            alc_final_flux = np.sum(norm_flux/norm_err**2,axis=1) / np.sum(1/norm_err**2,axis=1)
            alc_final_err = np.sqrt( 1 / np.sum(1/norm_err**2,axis=1))    

            targ_flux_corr = targ_flux_norm/alc_final_flux
            targ_flux_corr_err = np.sqrt( (targ_flux_err_norm/alc_final_flux)**2 + (targ_flux_norm*alc_final_err/(alc_final_flux**2))**2 )

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
            
            #Plot the corrected target lightcurve. 
            pines_lc_plot(times, targ_flux_corr, binned_times, binned_flux, binned_errs, analysis_path, time_mode='UT')
            
            cmap = plt.cm.viridis
            fig, ax = plt.subplots(num_refs+1, 1, figsize=(14,5*num_refs))
            plt.subplots_adjust(top=0.95, left=0.05, bottom=0.05, right=0.95)

            ax[0].axhline(1, color='k', alpha=0.7, lw=1, zorder=0)
            ax[0].grid(alpha=0.4)
            ax[0].plot(times, targ_flux_corr, '.', zorder=1, color=cmap(0), alpha=0.5)
            ax[0].errorbar(binned_times, binned_flux, binned_errs, linestyle='', marker='o', zorder=2, color=cmap(0))
            ax[0].set_ylabel('Corrected Flux', fontsize=14)
            ax[0].set_title(short_name, fontsize=16, color=cmap(0))
            ax[0].tick_params(labelsize=12)
            ax[0].set_ylim(0.9,1.1)
            ax[0].set_xlabel('Time (JD$_{UTC})$', fontsize=14)
            
            for k in range(num_refs):
                color_ind = int(((k+2) * 256 / (num_refs+1))) - 1
                color = cmap(color_ind)
                ax[k+1].axhline(1, color='k', alpha=0.7, lw=1, zorder=0)
                ax[k+1].grid(alpha=0.4)
                ax[k+1].plot(times, corr_flux[:,k], '.', zorder=1, color=color, alpha=0.5)
                ax[k+1].errorbar(binned_times, binned_ref_fluxes[:,k], binned_ref_flux_errs[:,k], linestyle='', marker='o', zorder=2, color=color)
                ax[k+1].set_ylabel('Corrected Flux', fontsize=14)
                ax[k+1].set_title('CS '+str(k+1), fontsize=16, color=color)
                ax[k+1].tick_params(labelsize=12)
                ax[k+1].set_ylim(0.9,1.1)
                ax[k+1].set_xlabel('Time (JD$_{UTC})$', fontsize=14)
            plt.tight_layout()
            plt.savefig(analysis_path/'cs_fig.png')
            
            bin_dict = {'Binned Time':binned_times, 'Binned '+short_name+ ' Flux':binned_flux, 'Binned '+short_name+ ' Flux Error':binned_errs}
            for k in range(num_refs):
                bin_dict['Binned '+ref_names[k]+ ' Flux'] = binned_ref_fluxes[:,k]
                bin_dict['Binned '+ref_names[k]+ ' Flux Error'] = binned_ref_flux_errs[:,k]

            targ_dict = {'Time':times, 'Flux':targ_flux_corr}
