import numpy as np
import matplotlib.pyplot as plt
import pickle
import pdb
from scipy.stats.stats import pearsonr
from sklearn import linear_model
from astropy.stats import sigma_clipped_stats
from pandas import DataFrame
import copy 
import itertools
import julian
import matplotlib
from pylab import *
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.analysis.night_splitter import night_splitter
from pines_analysis_toolkit.analysis.block_splitter import block_splitter
import natsort
import pandas as pd
from scipy.stats import sigmaclip
import matplotlib.dates as mdates
from photutils import make_source_mask 

#Input parameters
def lightcurve(target, sources, centroided_sources, phot_type='aper', ref_set_choice=[], plot_mode='combined'):
    '''Authors: 
        Patrick Tamburo, Boston University, June 2020
    Purpose: 
            Makes a "simple" lightcurve with each reference star weighted equally when creating the artificial comparison lightcurve. 
    Inputs:
        target (str): the target's long name (e.g. '2MASS J12345678+1234567').
        sources (pandas DataFrame): DataFrame with source names, and x/y positions in the source_detect_image. Output from ref_star_chooser. 
        centroided_sources (pandas DataFrame): Dataframe with source positions in every image. Output from centroider. 
        phot_type (str): 'aper' or 'psf'. Whether to use aperture or PSF phomometry NOTE: PSF photometry currently not implemented.
        ref_set_choice (list): list of reference IDs to use to make the lightcurve, in case you want to exclude any. 
        plot_mode (str): 'combined' or 'separate'. 'combined' plots all nights in one figure, while 'separate' plots nights in separate figures. 
    Outputs:
        Saves lightcurve plots to target's analysis directory. 
    TODO:
        PSF photometry 
        Regression? 
'''

    def regression(flux, seeing, airmass, corr_significance=1e-5):
        #Looks at correlations between seeing and airmass with the target flux.
        #Takes those variables which are significantly correlated, and uses them in a linear regression to de-correlated the target flux. 
        
        #Use the seeing in the regression if it's significantly correlated
        if pearsonr(seeing,flux)[1] < corr_significance:
            use_seeing = True
        else:
            use_seeing = False
        
        #Same thing for airmass
        if pearsonr(airmass,flux)[1] < corr_significance:
            use_airmass = True
        else:
            use_airmass = False   
        
        #Now set up the linear regression.
        regr = linear_model.LinearRegression()
        
        regress_dict = {}

        #Add seeing, background, and airmass, if applicable.
        if use_seeing:
            key = 'seeing'
            regress_dict[key] = seeing
        
        if use_airmass:
            key = 'airmass'
            regress_dict[key] = airmass
        
        #Finally, add target flux
        regress_dict['flux'] = flux

        #Get list of keys
        keylist = list()
        for i in regress_dict.keys():
            keylist.append(i)


        #Create data frame of regressors.
        df = DataFrame(regress_dict,columns=keylist)
        x = df[keylist[0:len(keylist)-1]]
        y = df['flux']

        if np.shape(x)[1] >0:
            regr.fit(x,y)

            #Now, define the model.
            linear_regression_model = regr.intercept_
            i = 0

            #Add in the other regressors, as necessary. Can't think of a way of doing this generally, just use a bunch of ifs. 
            if (use_seeing) and (use_airmass):
                linear_regression_model = linear_regression_model + regr.coef_[0]*seeing + regr.coef_[1]*airmass
            if (use_seeing) and not (use_airmass):
                linear_regression_model = linear_regression_model + regr.coef_[0]*seeing
            if not (use_seeing) and (use_airmass):
                linear_regression_model = linear_regression_model + regr.coef_[0]*airmass
            #Divide out the fit. 
            corrected_flux = flux/linear_regression_model    
        else:
            #print('No regressors used.')
            corrected_flux = flux
        return corrected_flux
    plt.ion()
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    outlier_tolerance = 0.2 #If a reference > outlier_tolerance of its values above sigma clipping threshold, mark it as bad. 
    centroided_sources.columns = centroided_sources.keys().str.strip()
    
    #Get list of photometry files for this target. 
    photometry_path = pines_path/('Objects/'+short_name+'/'+phot_type+'_phot/')
    analysis_path = pines_path/('Objects/'+short_name+'/analysis')
    photometry_files = natsort.natsorted([x for x in photometry_path.glob('*.csv')])

    num_refs = len(sources) - 1
    
    #Loop over all photometry files in the aper_phot directory. 
    for i in range(len(photometry_files)):
        #Load in the photometry data. 
        if phot_type =='aper':
            aperture_radius = float(str(photometry_files[i]).split('_')[-3])
        
        phot_data = pines_log_reader(photometry_files[i])

        #Remove entries that have NaN's for flux values.
        for j in range(len(sources['Name'])):
            name = sources['Name'][j]
            phot_data[name+' Flux'] = phot_data[name+' Flux'].astype(float)
            phot_data[name+' Flux Error'] = phot_data[name+' Flux Error'].astype(float)

        #Get target interpolation warnings. 
        targ_interp_flags = np.array(phot_data[short_name+' Interpolation Flag'])

        #Get times of exposures. 
        times = np.array(phot_data['Time JD'])
        seeing = np.array(phot_data['Seeing'])
        airmass = np.array(phot_data['Airmass'])
        background = np.array(phot_data[sources['Name'][0]+' Background'])
        
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
            values, clow, chigh = sigmaclip(ref_flux[j], low=2.5, high=2.5)
            if (len(phot_data) - len(values)) > (int(outlier_tolerance * len(phot_data))):
                print('Have to add flagging bad refs.')

        closest_ref = np.where(abs(np.nanmean(ref_flux, axis=1)-np.nanmean(targ_flux)) == min(abs(np.nanmean(ref_flux, axis=1)-np.nanmean(targ_flux))))[0][0]

        #Split data up into individual nights. 
        night_inds = night_splitter(times)
        num_nights = len(night_inds)
        
        if plot_mode == 'combined':
            fig, axis = plt.subplots(nrows=1, ncols=num_nights, figsize=(16, 5))
        
        colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

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

        for j in range(num_nights):
            if plot_mode == 'separate':
                fig, ax = plt.subplots(1, 1, figsize=(16,5))
            else:
                ax = axis[j]
            
            if phot_type =='aper':
                fig.suptitle(short_name, fontsize=16)
           
            j += 1
            if j == 1:
                ax.set_ylabel('Normalized Flux', fontsize=16)
            ax.set_xlabel('Time (UT)', fontsize=16)

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

            #Correct the example reference lightcurve using the alc. 
            ref_corr_norm = ref_flux[closest_ref][inds] / np.nanmedian(ref_flux[closest_ref][inds])
            ref_corr = ref_corr_norm / alc
            ref_corr = ref_corr / np.nanmedian(ref_corr)

            #Plot the target and reference lightcurves. 
            t_plot, = ax.plot(dts[inds], targ_corr, '.', color=colors[i])
            line_list.append(t_plot)
            myFmt = mdates.DateFormatter('%H:%M')
            ax.xaxis.set_major_formatter(myFmt)
            fig.autofmt_xdate()       

            #Do sigma clipping on the corrected lightcurve to get rid of outliers (from clouds, bad target centroid, cosmic rays, etc.)
            ###vals, lo, hi = sigmaclip(targ_corr, low=2.5, high=2.5)
            avg, med, std = sigma_clipped_stats(targ_corr, sigma=3)      
            bad_vals = np.where((targ_corr > med + 5*std) | (targ_corr < med - 5*std))[0]
            good_vals = np.where((targ_corr < med + 5*std) & (targ_corr > med - 5*std))[0]      
            vals = targ_corr[good_vals]
            if len(bad_vals) != 0:
                plt.plot(dts[inds][bad_vals], targ_corr[bad_vals], marker='x',color='r', mew=1.8, ms=7, zorder=0, ls='')
            
            blocks = block_splitter(times[inds], bad_vals)
            bin_times = np.zeros(len(blocks))
            bin_fluxes = np.zeros(len(blocks))
            bin_errs = np.zeros(len(blocks))
            bin_dts = []
            for k in range(len(blocks)):
                try: 
                    bin_times[k] = np.mean(times[inds][blocks[k]])
                    #vals, hi, lo = sigmaclip(targ_corr[blocks[k]],high=3,low=3) #Exclude outliers. 
                    bin_fluxes[k] = np.mean(targ_corr[blocks[k]])
                    bin_errs[k] = np.std(targ_corr[blocks[k]]) / np.sqrt(len(targ_corr[blocks[k]]))
                    bin_dts.append(julian.from_jd(bin_times[k], fmt='jd'))
                except:
                    pdb.set_trace()
            bin_dts = np.array(bin_dts)
            ax.errorbar(bin_dts, bin_fluxes, yerr=bin_errs, marker='o', color='k',zorder=3, ls='')
            
            #Draw the y=1 and 5-sigma detection threshold lines. 
            ax.axhline(y=1, color='r', lw=2, zorder=0)
            ax.axhline(1-5*np.median(bin_errs), zorder=0, lw=2, color='k', ls='--', alpha=0.4)

            #Set the y-range so you can see the 5-sigma detection line. 
            ax.set_ylim(0.9, 1.1)

            #Set the x-range to be the same for all nights. 
            ax.set_xlim(julian.from_jd(times[inds][0]-0.025, fmt='jd'), julian.from_jd(times[inds][0]+longest_night+0.025, fmt='jd'))
            
            ax.grid(alpha=0.2)
            ax.set_title(phot_data['Time UT'][inds[0]].split('T')[0], fontsize=14)
            ax.tick_params(labelsize=12)
            print('average seeing, night {}: {}'.format(j, np.mean(seeing[inds])))
            pdb.set_trace()
            #print('pearson correlation between target and closest ref: {}'.format(pearsonr(targ_corr[good_vals], ref_corr[good_vals])))
            
        print(np.mean(bin_errs))
        print('')
        fig.tight_layout(rect=[0, 0.03, 1, 0.93])
        if phot_type == 'aper':
            plt.savefig(analysis_path/(short_name+'_simple_lc_'+str(np.round(aperture_radius,1))+'.png'))
        
        pdb.set_trace()