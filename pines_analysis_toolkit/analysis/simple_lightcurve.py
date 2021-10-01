import numpy as np
import matplotlib.pyplot as plt
import pdb
from astropy.stats import sigma_clipped_stats
import julian
from pylab import *
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader
from pines_analysis_toolkit.analysis.night_splitter import night_splitter
from pines_analysis_toolkit.analysis.block_splitter import block_splitter
from pines_analysis_toolkit.utils import get_source_names
import natsort
import pandas as pd
import os.path

#Input parameters
def simple_lightcurve(short_name, phot_type='aper'):
    """Makes a 'simple' lightcurve with each reference star weighted equally when creating the artificial comparison lightcurve. Makes a light curve csv for every photometry file of type phot_type.

    :param short_name: short name of the target
    :type target: str
    :param phot_type: type of photometry, 'aper' or 'psf', defaults to 'aper'
    :type phot_type: str, optional
    """
    print('\nRunning simple_lightcurve().\n')

    pines_path = pines_dir_check()

    outlier_tolerance = 0.2 #If a reference > outlier_tolerance of its values above sigma clipping threshold, mark it as bad. 
    #centroided_sources.columns = centroided_sources.keys().str.strip()

    centroid_path = pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv')
    if not os.path.exists(centroid_path):
        raise RuntimeError('Centroid csv does not exist for {}, you need to make one using centroider.'.format(short_name))
    centroided_sources = pines_log_reader(centroid_path)
    sources = get_source_names(centroided_sources)

    #Get list of photometry files for this target. 
    photometry_path = pines_path/('Objects/'+short_name+'/'+phot_type+'_phot/')
    analysis_path = pines_path/('Objects/'+short_name+'/analysis/simple_lc_analysis/')
    if not os.path.exists(analysis_path):
        os.mkdir(analysis_path)
    photometry_files = natsort.natsorted([x for x in photometry_path.glob('*.csv')])

    num_refs = len(sources) - 1
    
    #Loop over all photometry files in the aper_phot directory. 
    for i in range(len(photometry_files)):
        #Load in the photometry data. 
        if phot_type =='aper':
            aperture_radius = float(str(photometry_files[i]).split('_')[-3])
        
        phot_data = pines_log_reader(photometry_files[i])

        #Remove entries that have NaN's for flux values.
        for j in range(len(sources)):
            name = sources[j]
            phot_data[name+' Flux'] = phot_data[name+' Flux'].astype(float)
            phot_data[name+' Flux Error'] = phot_data[name+' Flux Error'].astype(float)

        #Get target interpolation warnings. 
        targ_interp_flags = np.array(phot_data[short_name+' Interpolation Flag'])

        #Get times of exposures. 
        times = np.array(phot_data['Time BJD TDB'])
        
        #Convert to datetimes for plotting purposes.
        dts = np.array([julian.from_jd(times[i], fmt='jd') for i in range(len(times))])
        
        #Get the target's flux and background
        targ_flux = np.array(phot_data[short_name+' Flux'])
        targ_flux_err = np.array(phot_data[short_name+' Flux Error'])

        #Get the reference stars' fluxes and backgrounds. 
        ref_flux = np.zeros((num_refs, len(phot_data)))
        ref_flux_err = np.zeros((num_refs, len(phot_data)))
        for j in range(0, num_refs):
            ref_flux[j,:] = phot_data['Reference '+str(j+1)+' Flux']
            ref_flux_err[j,:] = phot_data['Reference '+str(j+1)+' Flux Error']
            #Discard variable stars. 
            #values, clow, chigh = sigmaclip(ref_flux[j], low=2.5, high=2.5)
            # if (len(phot_data) - len(values)) > (int(outlier_tolerance * len(phot_data))):
            #     print('Have to add flagging bad refs.')

        closest_ref = np.where(abs(np.nanmean(ref_flux, axis=1)-np.nanmean(targ_flux)) == min(abs(np.nanmean(ref_flux, axis=1)-np.nanmean(targ_flux))))[0][0]

        #Split data up into individual nights. 
        night_inds = night_splitter(times)
        num_nights = len(night_inds)
        
        #if plot_mode == 'combined':
        #    fig, axis = plt.subplots(nrows=1, ncols=num_nights, figsize=(16, 5))
        
        #colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

        #Get the time range of each night. Set each plot panel's xrange according to the night with the longest time. Makes seeing potential variability signals easier. 
        night_lengths = np.zeros(num_nights)
        for j in range(num_nights):
            inds = night_inds[j]
            night_lengths[j] = times[inds][-1] - times[inds][0]
        longest_night = max(night_lengths)
        longest_night_hours = np.ceil(longest_night*24)
        global line_list, filename_list
        line_list = []
        filename_list = []
        flux_save = []
        for j in range(num_nights):
            j += 1

            inds = night_inds[j-1]
            filename_list.append(np.array([phot_data['Filename'][z] for z in inds]))
            alc = np.zeros(len(inds))

            #Normalize reference lightcurves
            #TODO: Each night should be normalized separately. 
            for k in range(num_refs):
                ref_flux[k][inds] = ref_flux[k][inds] / np.nanmedian(ref_flux[k][inds])

            for k in range(len(inds)):
                #Do a sigma clip on normalized references to avoid biasing median. 
                ###values, clow, chigh = sigmaclip(ref_flux[:,inds[k]][~np.isnan(ref_flux[:,inds[k]])], low=1.5, high=1.5) 
                ###alc[k] = np.median(values)
                avg, med, std = sigma_clipped_stats(ref_flux[:,inds[k]][~np.isnan(ref_flux[:,inds[k]])], sigma=1.5)
                alc[k] = med
            
            #Correct the target lightcurve using the alc. 
            alc = alc / np.nanmedian(alc)
            targ_flux_norm = targ_flux[inds] / np.nanmedian(targ_flux[inds])
            targ_corr = targ_flux_norm / alc
            targ_corr = targ_corr / np.nanmedian(targ_corr)
            flux_save.extend(targ_corr)

            #Correct the example reference lightcurve using the alc. 
            ref_corr_norm = ref_flux[closest_ref][inds] / np.nanmedian(ref_flux[closest_ref][inds])
            ref_corr = ref_corr_norm / alc
            ref_corr = ref_corr / np.nanmedian(ref_corr)

            # #Do sigma clipping on the corrected lightcurve to get rid of outliers (from clouds, bad target centroid, cosmic rays, etc.)
            # ###vals, lo, hi = sigmaclip(targ_corr, low=2.5, high=2.5)
            # avg, med, std = sigma_clipped_stats(targ_corr, sigma=3)      
            # bad_vals = np.where((targ_corr > med + 5*std) | (targ_corr < med - 5*std))[0]
            # good_vals = np.where((targ_corr < med + 5*std) & (targ_corr > med - 5*std))[0]      

        #Output the simple lc data to a csv. 
        time_save = times        
        flux_err_save = np.zeros(len(flux_save)) + np.nanstd(targ_corr)
        output_dict = {'Time':time_save, 'Flux':flux_save, 'Flux Error':flux_err_save}
        output_df = pd.DataFrame(data=output_dict)
        if phot_type == 'aper':
            output_filename = analysis_path/(short_name+'_simple_lc_aper_phot_'+str(np.round(aperture_radius,1))+'_pix.csv')
            print('\nSaving to {}.\n'.format(output_filename))
            output_df.to_csv(output_filename)
        elif phot_type == 'psf':
            print("ERROR: Need to create flux output for PSF photometry.")     
               
    return