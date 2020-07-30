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

#Input parameters
def lightcurve(target, sources, phot_type='aper', ref_set_choice=[]):

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

    def noise_estimator(r,ref_set):
        G = 8.21 #e- / DN
        #RN = 19 #e- / pix. From Clemens (2007).
        RN = 32.89600862127023 #e- / pix. Measured on a set of biases from Nov. 2019.
        t = 30 #s
        npix = np.pi*r**2 #pix
        D = 0.98 #e- / pix / s

        ref_variances = np.zeros((len(ref_set),len(times))) #Record each reference star's contribution to the error. 
        ref_rates = np.zeros((len(ref_set),len(times)))
        ref_backgrounds = np.zeros((len(ref_set),len(times)))

        R_star = targ_flux * G / t #e- / s for the target. 
        R_ref  = ref_flux  * G / t #e- / s for the reference set. 

        for ii in range(len(ref_set)):
            ref_rates[ii] = restored_output['Flux'][ref_set][ii][use_locs_2]*G / t
            ref_backgrounds[ii] = restored_output['Background'][ref_set][ii][use_locs_2]*G/t
            ref_variances[ii] = ref_rates[ii]*t + ref_backgrounds[ii]*t*npix + RN**2*npix + D*npix*t
        
        total_ref_variance = np.sum(ref_variances,axis=0)

    plt.ion()
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    outlier_tolerance = 0.2 #If a reference > outlier_tolerance of its values above sigma clipping threshold, mark it as bad. 

    #Get list of photometry files for this target. 
    photometry_path = pines_path/('Objects/'+short_name+'/'+phot_type+'_phot/')
    analysis_path = pines_path/('Objects/'+short_name+'/analysis')
    photometry_files = natsort.natsorted([x for x in photometry_path.glob('*.csv')])

    num_refs = len(sources) - 1

    #Loop over all photometry files in the aper_phot directory. 
    for i in range(len(photometry_files)):
        #Load in the photometry data. 
        aperture_radius = float(str(photometry_files[i]).split('_')[-3])
        phot_data = pines_log_reader(photometry_files[i])

        #Get times of exposures. 
        times = np.array(phot_data['Time JD'])
        seeing = np.array(phot_data['Seeing'])
        airmass = np.array(phot_data['Airmass'])
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
                #pdb.set_trace()
        
        closest_ref = np.where(abs(np.mean(ref_flux, axis=1)-np.mean(targ_flux)) == min(abs(np.mean(ref_flux, axis=1)-np.mean(targ_flux))))[0][0]

        #Normalize reference lightcurves
        for k in range(num_refs):
            ref_flux[k] = ref_flux[k] / np.median(ref_flux[k])

        #Split data up into individual nights. 
        night_inds = night_splitter(times)
        num_nights = len(night_inds)
        fig = plt.figure(figsize=(16, 5))
        fig.suptitle(short_name+', aperture radius = '+str(np.round(aperture_radius,1))+' pixels', fontsize=16)
        colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

        #Get the time range of each night. Set each plot panel's xrange according to the night with the longest time. Makes seeing potential variability signals easier. 
        night_lengths = np.zeros(num_nights)
        for j in range(num_nights):
            inds = night_inds[j]
            night_lengths[j] = times[inds][-1] - times[inds][0]
        longest_night = max(night_lengths)
        longest_night_hours = np.ceil(longest_night*24)
        global ax_list, line_list, filename_list
        ax_list = []
        line_list = []
        filename_list = []
        for j in range(num_nights):
            j += 1
            ax = subplot(1, num_nights, j)
            ax_list.append(ax)
            if j == 1:
                ax.set_ylabel('Normalized Flux', fontsize=16)
            #else:
            #    ax.get_yaxis().set_visible(False)
            ax.set_xlabel('Time (UT)', fontsize=16)

            inds = night_inds[j-1]
            filename_list.append(np.array([phot_data['Filename'][z] for z in inds]))
            alc = np.zeros(len(inds))
            # if j == 5:
            #        plt.figure()
            for k in range(len(inds)):
                #Do a sigma clip on normalized references to avoid biasing median. 
                ###values, clow, chigh = sigmaclip(ref_flux[:,inds[k]][~np.isnan(ref_flux[:,inds[k]])], low=1.5, high=1.5) 
                ###alc[k] = np.median(values)
                avg, med, std = sigma_clipped_stats(ref_flux[:,inds[k]][~np.isnan(ref_flux[:,inds[k]])], sigma=1.5)
                alc[k] = med

                # if j == 5:
                #     plt.plot(np.zeros(len(ref_flux[:,inds[k]])) + times[inds[k]], ref_flux[:,inds[k]], 'k.')
                #     plt.plot(np.zeros(len(values)) + times[inds[k]], values, '.')
                #     plt.plot(times[inds[k]], alc[k], color='r', marker='o')
                #     pdb.set_trace()

            # if j == 5:
            #     pdb.set_trace()
            
            #Correct the target lightcurve using the alc. 
            alc = alc / np.median(alc)
            targ_flux_norm = targ_flux[inds] / np.median(targ_flux[inds])
            targ_corr = targ_flux_norm / alc
            targ_corr = targ_corr / np.median(targ_corr)

            
          
            #Correct the example reference lightcurve using the alc. 
            ref_corr_norm = ref_flux[closest_ref][inds] / np.median(ref_flux[closest_ref][inds])
            ref_corr = ref_corr_norm / alc
            ref_corr = ref_corr / np.median(ref_corr)


            #if j == 5:
            #    plt.plot(times[inds], targ_flux_norm, color='b', marker='o')
            #    pdb.set_trace()

            #Plot the target and reference lightcurves. 
            ax.plot(dts[inds], ref_corr,'.', color='tab:grey')
            t_plot, = ax.plot(dts[inds], targ_corr, '.', color=colors[i])
            line_list.append(t_plot)
            myFmt = mdates.DateFormatter('%H:%M')
            ax.xaxis.set_major_formatter(myFmt)
            fig.autofmt_xdate()       

            #Do sigma clipping on the corrected lightcurve to get rid of outliers (from clouds, bad target centroid, cosmic rays, etc.)
            ###vals, lo, hi = sigmaclip(targ_corr, low=2.5, high=2.5)
            avg, med, std = sigma_clipped_stats(targ_corr, sigma=3)      
            bad_vals = np.where((targ_corr > med + 3*std) | (targ_corr < med - 3*std))[0]
            good_vals = np.where((targ_corr < med + 3*std) & (targ_corr > med - 3*std))[0]      
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
            #print('pearson correlation between target and closest ref: {}'.format(pearsonr(targ_corr[good_vals], ref_corr[good_vals])))
            
        print(np.mean(bin_errs))
        print('')
        fig.tight_layout(rect=[0, 0.03, 1, 0.93])
        plt.savefig(analysis_path/(short_name+'_simple_lc_'+str(np.round(aperture_radius,1))+'.png'))
        pdb.set_trace()
        # #Weight reference stars, following Murray et al. (2020): https://ui.adsabs.harvard.edu/abs/2020MNRAS.495.2446M/abstract
        # #Weights based on distance from target. 
        # targ_loc = [sources['Source Detect X'][0], sources['Source Detect Y'][0]]
        # distances = np.array((sources[['Source Detect X', 'Source Detect Y']] - np.array(targ_loc)).pow(2).sum(1).pow(0.5))[1:] #Distances between references and targets, in pixels. 
        # s_max = max(distances)
        # a = 1 #A "parameter optimized each night" to minimize the average spread of the target's lightcurve... not sure what to use. 
        # W_dist = 1/(1+(a*distances/s_max)**2)
        # W_dist = W_dist / W_dist.sum() #Normalize to sum to 1.

        # #Loop over each night
        # for j in range(num_nights):
        #     inds = night_inds[j]
        #     #Weights based on photon noise.
        #     W_var = 1/(np.mean(ref_flux_err[:,inds], axis=1))**2
        #     W_var = W_var / W_var.sum() #Normalize to sum to 1. 

        #     max_difference = 9999. #Initialize. 

        #     while max_difference > 0.0000001:
        #         #Get total weights
        #         W = W_var * W_dist #The total weight given to each reference is just the product of its photon noise weight and its distance weight. 
        #         W = W / W.sum() #Normalize to sum to 1. 
        #         print(W)
        #         #time.sleep(2)
        #         print('')
        #         #Make the artifical light curve (alc).
        #         alc = np.zeros(len(inds))
        #         for k in range(len(inds)):
        #             alc[k] = sum(W * ref_flux[:,inds[k]]) / sum(W)
                
        #         #for k in range(num_refs):
        #         #   plt.plot(times[inds], ref_flux[k, inds]/np.median(ref_flux[k,inds]),marker='.')
        #         #plt.plot(times[inds], alc/np.median(alc), marker='o',color='k')
        #         #plt.plot(times[inds], targ_flux[inds]/np.median(targ_flux[inds]), marker='o',color='r')
        #         #pdb.set_trace()
        #         #plt.pause(0.5)
        #         #plt.clf()

        #         #REF FLUX NEEDS TO BE ABSOLUTE
        #         #Divide every reference's absolute lightcurve by this alc to produce a differential lightcurve. 
        #         ref_flux_diff = np.zeros((num_refs, len(inds)))
        #         W_var_new = np.zeros(num_refs)
                
        #         for k in range(num_refs):
        #             ref_flux_diff[k] = ref_flux[k,inds] / alc
        #             W_var_new[k] = 1/(np.std(ref_flux_diff[k])**2)
                   
        #         W_var_new = W_var_new / W_var_new.sum()
                
        #         max_difference = max(abs(W_var - W_var_new))
        #         print(max_difference)
        #         W_var = W_var_new
        #     plt.plot(times[inds],(targ_flux[inds]/alc)/np.median(targ_flux[inds]/alc),'.')
            
            # plt.figure()
            # for k in range(num_refs):
            #     plt.plot(times[inds], ref_flux[k, inds]/np.median(ref_flux[k,inds]),marker='.')
            # plt.plot(times[inds], alc/np.median(alc), marker='o',color='k')
            # plt.plot(times[inds], targ_flux[inds]/np.median(targ_flux[inds]), marker='o',color='r')
            #pdb.set_trace()
        

    # use_refs = '' #'all' for all refs, '' to optimize
    # bin_lc = 1
    # plots = 1
    # corr_significance = 1e-3 #For use in choosing which variables to use in linear regression. 
    # refs_to_remove = []
    # weight_ref_stars = 'y'
    # if (len(ref_set_choice) > 0):
    #     use_refs = 'choice'
    #     choice_set = [2,4,5,7,8,10]        
    #     choice_set = ref_set_choice


    # #Set these to zero to use all blocks in error analysis. 
    # #block_to_ignore_start = 2458883.7209662
    # #block_to_ignore_end = 2458883.73
    # block_to_ignore_start = 0
    # block_to_ignore_end = 0


    # #Pathing
    # local_path = '/Users/tamburo/Documents/Data/PINES/'
    # path = local_path+'Objects/'+target_name+'/aper_phot/photometry/'+target_name+'_ap='+aperture_radius+'_photometry.p'

    # #Read in output. 
    # restored_output = pickle.load(open(path,'rb'))

    # #Remove bad reference stars
    # restored_output['Flux'] = np.delete(restored_output['Flux'], refs_to_remove, axis=0)
    # restored_output['X position'] = np.delete(restored_output['X position'], refs_to_remove, axis=0)
    # restored_output['Y position'] = np.delete(restored_output['Y position'], refs_to_remove, axis=0)
    # restored_output['X subframe position'] = np.delete(restored_output['X subframe position'], refs_to_remove, axis=0)
    # restored_output['Y subframe position'] = np.delete(restored_output['Y subframe position'], refs_to_remove, axis=0)

    # #Throw out flagged images. 
    # times = restored_output['Time']

    # targ_flux = restored_output['Flux'][0,:]
    # x_position = restored_output['X position']
    # y_position = restored_output['Y position']
    # x_subframe_position = restored_output['X subframe position']
    # y_subframe_position = restored_output['Y subframe position']
    # targ_background = restored_output['Background'][0]
    # x_seeing = restored_output['X seeing']
    # y_seeing = restored_output['Y seeing']
    # airmass = restored_output['Airmass']

    # avg,med,std = sigma_clipped_stats(targ_flux)
    # #Make a cut on really faint images (this is my attempt to automatically flag cloudy images).
    # #use_locs_2 = np.where( (airmass < 2.5) & (targ_background < 6000) & (targ_background > 200) & (times > 2458885.3) & (targ_flux > (med-3*std)))[0] 
    # #use_locs_2 = np.where( (airmass < 2.5) & (targ_background < 3000) & (targ_background > 200) & (times > 2458883.3) & (times < 2458884.3) & (targ_flux > (med-3*std)))[0] 
    # use_locs_2 = np.where( (airmass < 2.5) & (targ_background < 3000) & (targ_background > 200)  & (targ_flux > (med-3*std)))[0] 
    # #use_locs_2 = np.where( (targ_flux/np.median(targ_flux) > 0.2) & (targ_background < 500))[0]

    # times = times[use_locs_2]
    # targ_flux = targ_flux[use_locs_2]
    # x_position = x_position[:,use_locs_2]
    # y_position = y_position[:,use_locs_2]
    # x_subframe_position = x_subframe_position[:,use_locs_2]
    # y_subframe_position = y_subframe_position[:,use_locs_2]
    # targ_background = targ_background[use_locs_2]
    # x_seeing = x_seeing[use_locs_2]
    # y_seeing = y_seeing[use_locs_2]
    # airmass = airmass[use_locs_2]

    # #Disregard a specified block of data from the standard deviation analysis, may bias your choice of reference stars.
    # if (block_to_ignore_start != 0) and (block_to_ignore_end != 0):
    #     norm_locs = np.where(( (times > times[0]) & (times < block_to_ignore_start))  | ((times > block_to_ignore_end) & (times < times[-1])))[0]
    # else:
    #     norm_locs = np.arange(0,len(times),1) #Else, use all the data. 
        
    # plt.ion()
    # fig,ax = plt.subplots(nrows=1,ncols=1,sharex=True,figsize=(7,6))
    # for i in range(len(restored_output['Flux'])):
    #     if i == 0:
    #         ax.plot(times,targ_flux,label='Target',lw=3,marker='o')
    #     else:
    #         ax.plot(times,restored_output['Flux'][i,:][use_locs_2],label='Ref. '+str(i),marker='o')
    # ax.legend()
    # plt.show()
    # # pdb.set_trace()

    # #Get all potential combinations of reference stars.
    # num_refs = np.shape(restored_output['Flux'][1:,:])[0] #Count the number of references
    # ref_numbers = np.arange(1,num_refs+1,1)
    # ref_combinations = []
    # for i in range(1,num_refs+1):
    #     sets = np.array(list(itertools.combinations(ref_numbers,i))).astype(int)
    #     for j in range(len(sets)):
    #         ref_combinations.append(sets[j])
    # ref_combinations = np.array(ref_combinations)
    # if use_refs == 'all': #Force all reference stars
    #     ref_combinations = [ref_combinations[-1],ref_combinations[-1]]

    # if use_refs == 'choice':
    #     ref_combinations = [choice_set,choice_set]

    # best_stdev = 100 #Initialize to large value.

    # #Loop over reference combinations to find the one with the smallest scatter. 
    # for i in range(len(ref_combinations)):
    #     ref_set = ref_combinations[i]

    #     if weight_ref_stars == 'n':
    #         ref_flux = []
            
    #         #First normalize the reference curves
    #         for j in range(len(ref_set)):
    #             ref_flux.append(restored_output['Flux'][ref_set][j][use_locs_2]/np.median(restored_output['Flux'][ref_set][j][use_locs_2]))
    #         ref_flux = np.mean(ref_flux,axis=0)
    #     else:
    #         ref_flux = sum(restored_output['Flux'][ref_set],axis=0)[use_locs_2]

    #     corrected_targ_flux = targ_flux/ref_flux

    #     norm_targ_flux = corrected_targ_flux/np.median(corrected_targ_flux)
    #     #Do the linear regression
    #     (norm_targ_flux_decorrelated,seeing) = regression()

    #     #print('Stddev. of lightcurve, post-regression: ',np.std(norm_targ_flux_decorrelated))
    #     std_plot =  np.std(norm_targ_flux_decorrelated[norm_locs])

    #     #Calculate the noise properties of the lightcurve: target, sky, reference, read noise, dark current. 
    #     noises = noise_estimator(float(aperture_radius),ref_set)

    #     if std_plot < best_stdev:
    #         print('New best set found: ',ref_set)
    #         best_set = ref_set
    #         best_targ_over_ref = np.mean(corrected_targ_flux)
    #         best_norm_targ_flux_decorrelated = norm_targ_flux_decorrelated
    #         best_ref_flux = ref_flux
    #         best_stdev = std_plot


    # print('')
    # print('Best set: ',best_set)
    # print('Best targ./ref.: ',best_targ_over_ref)


    # #Bin the data over blocks
    # if bin_lc:
    #     bin_locs = np.where(np.gradient(times)*24*60 > 10)[0]
    #     num_bins = int(len(bin_locs)/2+2) #Plus 2 from start/end bins
    #     binned_time = np.zeros(num_bins-1)
    #     binned_flux = np.zeros(num_bins-1)
    #     binned_err = np.zeros(num_bins-1)
    #     bin_starts = np.zeros(num_bins-1)
    #     bin_ends = np.zeros(num_bins-1)
    #     if len(bin_locs) > 0:
    #         for i in range(num_bins):
    #             if i == 0:
    #                 binned_time[i] = np.mean(times[0:bin_locs[0]])
    #                 binned_flux[i] = np.mean(best_norm_targ_flux_decorrelated[0:bin_locs[0]])
    #                 #binned_err[i] = np.std(best_norm_targ_flux_decorrelated[0:bin_locs[0]])/np.sqrt(bin_locs[0])
    #                 binned_err[i] = np.std(best_norm_targ_flux_decorrelated)/np.sqrt(bin_locs[0])
    #                 bin_starts[i] = 0
    #                 bin_ends[i] = bin_locs[0]
    #             elif i < num_bins-2:
    #                 binned_time[i] = np.mean(times[bin_locs[2*(i-1)+1]:bin_locs[2*(i-1)+2]])
    #                 binned_flux[i] = np.mean(best_norm_targ_flux_decorrelated[bin_locs[2*(i-1)+1]:bin_locs[2*(i-1)+2]])
    #                 #binned_err[i] = np.std(best_norm_targ_flux_decorrelated[bin_locs[2*(i-1)+1]:bin_locs[2*(i-1)+2]])/np.sqrt(bin_locs[2*(i-1)+2] - bin_locs[2*(i-1)+1])
    #                 binned_err[i] = np.std(best_norm_targ_flux_decorrelated)/np.sqrt(bin_locs[2*(i-1)+2] - bin_locs[2*(i-1)+1])
    #                 bin_starts[i] = bin_locs[2*(i-1)+1]
    #                 bin_ends[i] = bin_locs[2*(i-1)+2]
    #             elif i == num_bins-2:
    #                 binned_time[i] = np.mean(times[bin_locs[-1]:len(times)])
    #                 binned_flux[i] = np.mean(best_norm_targ_flux_decorrelated[bin_locs[-1]:len(times)])
    #                 #binned_err[i] = np.std(best_norm_targ_flux_decorrelated[bin_locs[-1]:len(times)])/np.sqrt(len(times) - bin_locs[-1])
    #                 binned_err[i] = np.std(best_norm_targ_flux_decorrelated)/np.sqrt(len(times) - bin_locs[-1])
    #                 bin_starts[i] = bin_locs[-1]
    #                 bin_ends[i] = len(times)

    #     print('Best stdev: ',np.std(best_norm_targ_flux_decorrelated))
    #     print('Best mean binned err: ',np.mean(binned_err))


    # if plots:
    #     times_plot = [julian.from_jd(times[i], fmt='jd') for i in range(len(times))]
    #     if len(bin_locs ) > 0:
    #         binned_time_plot = [julian.from_jd(binned_time[i], fmt='jd') for i in range(len(binned_time))]

    #     #First just plot the lightcurve. 
    #     fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(10,4))
    #     ax.plot([min(times_plot),max(times_plot)],[1,1],color='r',lw=2,zorder=0)
    #     ax.plot(times_plot,best_norm_targ_flux_decorrelated,marker='o',color='#DCDCDC',linestyle='',zorder=1)
    #     if len(bin_locs) > 0:
    #         ax.errorbar(binned_time_plot,binned_flux,binned_err,linestyle='',marker='o',color='k',capsize=3)
    #     #ax.plot(binned_time_plot,binned_flux,linestyle='',marker='o',color='k')
    #     ax.plot([min(times_plot),max(times_plot)],[1-5*np.mean(binned_err),1-5*np.mean(binned_err)],color='k',linestyle='--',lw=2,alpha=0.7)
    #     ax.text(max(times_plot),1-5*binned_err[0]+.001,'5$\sigma$ threshold')
    #     ax.set_ylabel('Norm. Flux',fontsize=15)
    #     ax.set_xlabel('JD$_{UTC}$',fontsize=15)
    #     ax.grid(alpha=0.4)
    #     ax.get_yaxis().set_label_coords(-0.06,0.5)
    #     ax.set_xlim(min(times_plot),max(times_plot))
    #     plt.tight_layout()
    #     plt.ion()
    #     # thismanager = get_current_fig_manager()
    #     # thismanager.window.wm_geometry("-1500+900")
    #     plt.show()

    #     #Make diagnostics plot. 
    #     fig,ax = plt.subplots(nrows=7,ncols=1,figsize=(13,10),sharex=True,gridspec_kw={'height_ratios': [2.5,1,1,1,1,1,1]})
    #     ax[0].plot(times_plot,best_norm_targ_flux_decorrelated,'k.',zorder=1)
    #     ax[0].plot([min(times_plot),max(times_plot)],[1,1],color='r',lw=2,zorder=0)
    #     ax[0].set_ylabel('Norm. Flux',fontsize=9)
    #     ax[0].grid(alpha=0.4)
    #     ax[0].get_yaxis().set_label_coords(-0.06,0.5)

    #     ax[1].plot(times_plot,targ_flux/np.mean(targ_flux),label='Target',alpha=0.7,marker='.')
    #     ax[1].plot(times_plot,best_ref_flux/np.mean(best_ref_flux),label='Reference',alpha=0.7,marker='.')
    #     ax[1].legend(bbox_to_anchor=(0.595, 0.525, 0.5, 0.5),fontsize=7)
    #     ax[1].set_ylabel('Norm. Flux',fontsize=9)
    #     ax[1].get_yaxis().set_label_coords(-0.06,0.5)

    #     #if use_refs == '' or use_refs == 'choice':
    #     for i in range(len(best_set)+1):
    #         if i == 0:
    #             lab = 'Target'
    #         else:
    #             lab = 'Ref. '+str(best_set[i-1])
            
    #         ax[2].plot(times_plot,x_subframe_position[i],alpha=0.5,marker='.',label=lab)
    #         ax[3].plot(times_plot,y_subframe_position[i],alpha=0.5,marker='.')
    #     ax[2].set_ylabel('X sub. pos.',fontsize=9)
    #     ax[2].get_yaxis().set_label_coords(-0.06,0.5)
    #     ax[2].legend(bbox_to_anchor=(0.5775, 0.525, 0.5, 0.5),fontsize=7)

    #     ax[3].set_ylabel('Y sub. pos.',fontsize=9)
    #     ax[3].get_yaxis().set_label_coords(-0.06,0.5)


    #     ax[4].plot(times_plot,targ_background,marker='.')
    #     ax[4].set_ylabel('Bkg.',fontsize=9)
    #     ax[4].get_yaxis().set_label_coords(-0.06,0.5)

    #     ax[5].plot(times_plot,x_seeing,color='purple',marker='.',label='X seeing')
    #     ax[5].plot(times_plot,y_seeing,color='orange',marker='.',alpha=0.5,label='Y seeing')
    #     ax[5].legend()
    #     ax[5].set_ylabel('Seeing (")',fontsize=9)
    #     ax[5].set_xlim(min(times_plot),max(times_plot))
    #     ax[5].get_yaxis().set_label_coords(-0.06,0.5)

    #     ax[6].plot(times_plot,airmass,color='gold',marker='.',alpha=0.7)
    #     ax[6].set_ylabel('Airmass',fontsize=9)
    #     ax[6].set_xlabel('Julian Date',fontsize=9)
    #     ax[6].set_xlim(min(times_plot),max(times_plot))
    #     ax[6].get_yaxis().set_label_coords(-0.06,0.5)
    #     plt.suptitle(target_name)
    #     plt.ion()
    #     plt.subplots_adjust(top=0.95,left=0.08,bottom=0.05,right=0.92,hspace=0.1)
    #     plt.savefig(local_path+'Objects/'+target_name+'/analysis/'+target_name+'_diagnostics_plot.png')
        
    #     # #Make plot showing zoom-ins on all nights 
    #     # night_locs = np.where(np.gradient(times) > 0.2)[0]
    #     # num_nights = int(len(night_locs)/2)+1
    #     # fig,ax = plt.subplots(num_nights,1,figsize=(8,10))
    #     # for i in range(num_nights):
    #     #     if i == 0:
    #     #         ax[i].plot(times_plot[0:night_locs[0]],best_norm_targ_flux_decorrelated[0:night_locs[0]],'k.')
    #     #         ax[i].set_title(target_name,fontsize=16)
    #     #     elif i < num_nights - 1:
    #     #         ax[i].plot(times_plot[night_locs[2*(i-1)+1]:night_locs[2*(i-1)+2]],best_norm_targ_flux_decorrelated[night_locs[2*(i-1)+1]:night_locs[2*(i-1)+2]],'k.')
    #     #     else:
    #     #         ax[i].plot(times_plot[night_locs[-1]:],best_norm_targ_flux_decorrelated[night_locs[-1]:],'k.')
    #     # plt.tight_layout()

    #     save_dict = {'Times':times_plot,'Flux':best_norm_targ_flux_decorrelated,'Bin starts':bin_starts,'Bin ends':bin_ends}
    #     pickle.dump(save_dict,open(local_path+'Objects/'+target_name+'/analysis/'+target_name+'_lightcurve.p','wb'))
    #     plt.show()
    #     pdb.set_trace()


