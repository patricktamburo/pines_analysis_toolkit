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

#def lightcurve():

#Input parameters
def lightcurve(target_name,aperture_radius,ref_set_choice=[]):

    def regression():
        #Looks at correlations between seeing, background, and airmass with the target flux.
        #Takes those variables which are significantly correlated, and uses them in a linear regression to de-correlated the target flux. 
        

        
        #Calculate correlation with x and y seeing, and use the one with the higher correlation for regression. 

        x_seeing_corr = pearsonr(norm_targ_flux,x_seeing)[0]
        y_seeing_corr = pearsonr(norm_targ_flux,y_seeing)[0]
        if max(abs(x_seeing_corr),abs(y_seeing_corr)) == abs(x_seeing_corr):
            seeing = x_seeing
        elif max(abs(x_seeing_corr),abs(y_seeing_corr)) == abs(y_seeing_corr):
            seeing = y_seeing
        #TODO: have a bug where y_seeing is like 3 times larger than x_seeing, just use x_seeing for now.
        #seeing = x_seeing
        #Use the seeing in the regression if it's significantly correlated
        if pearsonr(seeing,norm_targ_flux)[1] < corr_significance:
            use_seeing = True
        else:
            use_seeing = False
        
        #Same thing for background
        if pearsonr(targ_background,norm_targ_flux)[1] < corr_significance:
            use_background = True
        else:
            use_background = False
        
        #Same thing for airmass
        if pearsonr(airmass,norm_targ_flux)[1] < corr_significance:
            use_airmass = True
        else:
            use_airmass = False   


        #print('Correlation with x seeing: ',x_seeing_corr)
        #print('Correlation with y seeing: ',y_seeing_corr)
        #print('Correlation with background: ',pearsonr(norm_targ_flux,targ_background)[0])
        #print('Correlation with airmass: ',pearsonr(norm_targ_flux,airmass)[0])
        
        #print('Stddev. of lightcurve, pre-regression: ',np.std(norm_targ_flux))

        #Now set up the linear regression.
        regr = linear_model.LinearRegression()
        
        regress_dict = {}

        #Add seeing, background, and airmass, if applicable.
        if use_seeing:
            key = 'seeing'
            regress_dict[key] = seeing
        
        if use_background:
            key = 'background'
            regress_dict[key] = targ_background
        
        if use_airmass:
            key = 'airmass'
            regress_dict[key] = airmass
        
        # if use_x_pos:
        #     key = 'x_pos'
        #     regress_dict[key] = x_position[0,:]

        #Finally, add target flux
        regress_dict['Norm targ flux'] = norm_targ_flux

        

        #Get list of keys
        keylist = list()
        for i in regress_dict.keys():
            keylist.append(i)


        #Create data frame of regressors.
        df = DataFrame(regress_dict,columns=keylist)
        x = df[keylist[0:len(keylist)-1]]
        y = df['Norm targ flux']

        if np.shape(x)[1] >0:
            regr.fit(x,y)

            #Now, define the model.
            linear_regression_model = regr.intercept_
            i = 0

            #Add in the other regressors, as necessary. Can't think of a way of doing this generally, just use a bunch of ifs. 
            
            if (use_seeing) and (use_background) and (use_airmass):
                linear_regression_model = linear_regression_model+regr.coef_[i]*seeing+regr.coef_[i+1]*targ_background+regr.coef_[i+2]*airmass
            if (use_seeing) and (use_background) and not (use_airmass):
                linear_regression_model = linear_regression_model+regr.coef_[i]*seeing+regr.coef_[i+1]*targ_background
            if (use_seeing) and (use_airmass) and not (use_background):
                linear_regression_model = linear_regression_model+regr.coef_[i]*seeing+regr.coef_[i+1]*airmass
            if (use_background) and (use_airmass) and not (use_seeing):
                linear_regression_model = linear_regression_model+regr.coef_[i]*targ_background+regr.coef_[i+1]*airmass
            if (use_seeing) and not (use_background) and not (use_airmass):
                linear_regression_model = linear_regression_model+regr.coef_[i]*seeing
            if (use_background) and not (use_seeing) and not (use_airmass):
                linear_regression_model = linear_regression_model+regr.coef_[i]*targ_background
            if (use_airmass) and not (use_seeing) and not (use_background):
                linear_regression_model = linear_regression_model+regr.coef_[i]*airmass
            

            
            #Divide out the fit. 
            corrected_flux = norm_targ_flux/linear_regression_model    
        else:
            #print('No regressors used.')
            corrected_flux = norm_targ_flux
        return (corrected_flux,seeing)

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


    use_refs = '' #'all' for all refs, '' to optimize
    bin_lc = 1
    plots = 1
    corr_significance = 1e-3 #For use in choosing which variables to use in linear regression. 
    refs_to_remove = []
    weight_ref_stars = 'y'
    if (len(ref_set_choice) > 0):
        use_refs = 'choice'
        choice_set = [2,4,5,7,8,10]        
        choice_set = ref_set_choice


    #Set these to zero to use all blocks in error analysis. 
    #block_to_ignore_start = 2458883.7209662
    #block_to_ignore_end = 2458883.73
    block_to_ignore_start = 0
    block_to_ignore_end = 0


    #Pathing
    local_path = '/Users/tamburo/Documents/Data/PINES/'
    path = local_path+'Objects/'+target_name+'/aper_phot/photometry/'+target_name+'_ap='+aperture_radius+'_photometry.p'

    #Read in output. 
    restored_output = pickle.load(open(path,'rb'))

    #Remove bad reference stars
    restored_output['Flux'] = np.delete(restored_output['Flux'], refs_to_remove, axis=0)
    restored_output['X position'] = np.delete(restored_output['X position'], refs_to_remove, axis=0)
    restored_output['Y position'] = np.delete(restored_output['Y position'], refs_to_remove, axis=0)
    restored_output['X subframe position'] = np.delete(restored_output['X subframe position'], refs_to_remove, axis=0)
    restored_output['Y subframe position'] = np.delete(restored_output['Y subframe position'], refs_to_remove, axis=0)

    #Throw out flagged images. 
    times = restored_output['Time']

    targ_flux = restored_output['Flux'][0,:]
    x_position = restored_output['X position']
    y_position = restored_output['Y position']
    x_subframe_position = restored_output['X subframe position']
    y_subframe_position = restored_output['Y subframe position']
    targ_background = restored_output['Background'][0]
    x_seeing = restored_output['X seeing']
    y_seeing = restored_output['Y seeing']
    airmass = restored_output['Airmass']

    avg,med,std = sigma_clipped_stats(targ_flux)
    #Make a cut on really faint images (this is my attempt to automatically flag cloudy images).
    #use_locs_2 = np.where( (airmass < 2.5) & (targ_background < 6000) & (targ_background > 200) & (times > 2458885.3) & (targ_flux > (med-3*std)))[0] 
    #use_locs_2 = np.where( (airmass < 2.5) & (targ_background < 3000) & (targ_background > 200) & (times > 2458883.3) & (times < 2458884.3) & (targ_flux > (med-3*std)))[0] 
    use_locs_2 = np.where( (airmass < 2.5) & (targ_background < 3000) & (targ_background > 200)  & (targ_flux > (med-3*std)))[0] 
    #use_locs_2 = np.where( (targ_flux/np.median(targ_flux) > 0.2) & (targ_background < 500))[0]

    times = times[use_locs_2]
    targ_flux = targ_flux[use_locs_2]
    x_position = x_position[:,use_locs_2]
    y_position = y_position[:,use_locs_2]
    x_subframe_position = x_subframe_position[:,use_locs_2]
    y_subframe_position = y_subframe_position[:,use_locs_2]
    targ_background = targ_background[use_locs_2]
    x_seeing = x_seeing[use_locs_2]
    y_seeing = y_seeing[use_locs_2]
    airmass = airmass[use_locs_2]

    #Disregard a specified block of data from the standard deviation analysis, may bias your choice of reference stars.
    if (block_to_ignore_start != 0) and (block_to_ignore_end != 0):
        norm_locs = np.where(( (times > times[0]) & (times < block_to_ignore_start))  | ((times > block_to_ignore_end) & (times < times[-1])))[0]
    else:
        norm_locs = np.arange(0,len(times),1) #Else, use all the data. 
        
    plt.ion()
    fig,ax = plt.subplots(nrows=1,ncols=1,sharex=True,figsize=(7,6))
    for i in range(len(restored_output['Flux'])):
        if i == 0:
            ax.plot(times,targ_flux,label='Target',lw=3,marker='o')
        else:
            ax.plot(times,restored_output['Flux'][i,:][use_locs_2],label='Ref. '+str(i),marker='o')
    ax.legend()
    plt.show()
    # pdb.set_trace()

    #Get all potential combinations of reference stars.
    num_refs = np.shape(restored_output['Flux'][1:,:])[0] #Count the number of references
    ref_numbers = np.arange(1,num_refs+1,1)
    ref_combinations = []
    for i in range(1,num_refs+1):
        sets = np.array(list(itertools.combinations(ref_numbers,i))).astype(int)
        for j in range(len(sets)):
            ref_combinations.append(sets[j])
    ref_combinations = np.array(ref_combinations)
    if use_refs == 'all': #Force all reference stars
        ref_combinations = [ref_combinations[-1],ref_combinations[-1]]

    if use_refs == 'choice':
        ref_combinations = [choice_set,choice_set]

    best_stdev = 100 #Initialize to large value.

    #Loop over reference combinations to find the one with the smallest scatter. 
    for i in range(len(ref_combinations)):
        ref_set = ref_combinations[i]

        if weight_ref_stars == 'n':
            ref_flux = []
            
            #First normalize the reference curves
            for j in range(len(ref_set)):
                ref_flux.append(restored_output['Flux'][ref_set][j][use_locs_2]/np.median(restored_output['Flux'][ref_set][j][use_locs_2]))
            ref_flux = np.mean(ref_flux,axis=0)
        else:
            ref_flux = sum(restored_output['Flux'][ref_set],axis=0)[use_locs_2]

        corrected_targ_flux = targ_flux/ref_flux

        norm_targ_flux = corrected_targ_flux/np.median(corrected_targ_flux)
        #Do the linear regression
        (norm_targ_flux_decorrelated,seeing) = regression()

        #print('Stddev. of lightcurve, post-regression: ',np.std(norm_targ_flux_decorrelated))
        std_plot =  np.std(norm_targ_flux_decorrelated[norm_locs])

        #Calculate the noise properties of the lightcurve: target, sky, reference, read noise, dark current. 
        noises = noise_estimator(float(aperture_radius),ref_set)

        if std_plot < best_stdev:
            print('New best set found: ',ref_set)
            best_set = ref_set
            best_targ_over_ref = np.mean(corrected_targ_flux)
            best_norm_targ_flux_decorrelated = norm_targ_flux_decorrelated
            best_ref_flux = ref_flux
            best_stdev = std_plot


    print('')
    print('Best set: ',best_set)
    print('Best targ./ref.: ',best_targ_over_ref)


    #Bin the data over blocks
    if bin_lc:
        bin_locs = np.where(np.gradient(times)*24*60 > 10)[0]
        num_bins = int(len(bin_locs)/2+2) #Plus 2 from start/end bins
        binned_time = np.zeros(num_bins-1)
        binned_flux = np.zeros(num_bins-1)
        binned_err = np.zeros(num_bins-1)
        bin_starts = np.zeros(num_bins-1)
        bin_ends = np.zeros(num_bins-1)
        if len(bin_locs) > 0:
            for i in range(num_bins):
                if i == 0:
                    binned_time[i] = np.mean(times[0:bin_locs[0]])
                    binned_flux[i] = np.mean(best_norm_targ_flux_decorrelated[0:bin_locs[0]])
                    #binned_err[i] = np.std(best_norm_targ_flux_decorrelated[0:bin_locs[0]])/np.sqrt(bin_locs[0])
                    binned_err[i] = np.std(best_norm_targ_flux_decorrelated)/np.sqrt(bin_locs[0])
                    bin_starts[i] = 0
                    bin_ends[i] = bin_locs[0]
                elif i < num_bins-2:
                    binned_time[i] = np.mean(times[bin_locs[2*(i-1)+1]:bin_locs[2*(i-1)+2]])
                    binned_flux[i] = np.mean(best_norm_targ_flux_decorrelated[bin_locs[2*(i-1)+1]:bin_locs[2*(i-1)+2]])
                    #binned_err[i] = np.std(best_norm_targ_flux_decorrelated[bin_locs[2*(i-1)+1]:bin_locs[2*(i-1)+2]])/np.sqrt(bin_locs[2*(i-1)+2] - bin_locs[2*(i-1)+1])
                    binned_err[i] = np.std(best_norm_targ_flux_decorrelated)/np.sqrt(bin_locs[2*(i-1)+2] - bin_locs[2*(i-1)+1])
                    bin_starts[i] = bin_locs[2*(i-1)+1]
                    bin_ends[i] = bin_locs[2*(i-1)+2]
                elif i == num_bins-2:
                    binned_time[i] = np.mean(times[bin_locs[-1]:len(times)])
                    binned_flux[i] = np.mean(best_norm_targ_flux_decorrelated[bin_locs[-1]:len(times)])
                    #binned_err[i] = np.std(best_norm_targ_flux_decorrelated[bin_locs[-1]:len(times)])/np.sqrt(len(times) - bin_locs[-1])
                    binned_err[i] = np.std(best_norm_targ_flux_decorrelated)/np.sqrt(len(times) - bin_locs[-1])
                    bin_starts[i] = bin_locs[-1]
                    bin_ends[i] = len(times)

        print('Best stdev: ',np.std(best_norm_targ_flux_decorrelated))
        print('Best mean binned err: ',np.mean(binned_err))


    if plots:
        times_plot = [julian.from_jd(times[i], fmt='jd') for i in range(len(times))]
        if len(bin_locs ) > 0:
            binned_time_plot = [julian.from_jd(binned_time[i], fmt='jd') for i in range(len(binned_time))]

        #First just plot the lightcurve. 
        fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(10,4))
        ax.plot([min(times_plot),max(times_plot)],[1,1],color='r',lw=2,zorder=0)
        ax.plot(times_plot,best_norm_targ_flux_decorrelated,marker='o',color='#DCDCDC',linestyle='',zorder=1)
        if len(bin_locs) > 0:
            ax.errorbar(binned_time_plot,binned_flux,binned_err,linestyle='',marker='o',color='k',capsize=3)
        #ax.plot(binned_time_plot,binned_flux,linestyle='',marker='o',color='k')
        ax.plot([min(times_plot),max(times_plot)],[1-5*np.mean(binned_err),1-5*np.mean(binned_err)],color='k',linestyle='--',lw=2,alpha=0.7)
        ax.text(max(times_plot),1-5*binned_err[0]+.001,'5$\sigma$ threshold')
        ax.set_ylabel('Norm. Flux',fontsize=15)
        ax.set_xlabel('JD$_{UTC}$',fontsize=15)
        ax.grid(alpha=0.4)
        ax.get_yaxis().set_label_coords(-0.06,0.5)
        ax.set_xlim(min(times_plot),max(times_plot))
        plt.tight_layout()
        plt.ion()
        # thismanager = get_current_fig_manager()
        # thismanager.window.wm_geometry("-1500+900")
        plt.show()

        #Make diagnostics plot. 
        fig,ax = plt.subplots(nrows=7,ncols=1,figsize=(13,10),sharex=True,gridspec_kw={'height_ratios': [2.5,1,1,1,1,1,1]})
        ax[0].plot(times_plot,best_norm_targ_flux_decorrelated,'k.',zorder=1)
        ax[0].plot([min(times_plot),max(times_plot)],[1,1],color='r',lw=2,zorder=0)
        ax[0].set_ylabel('Norm. Flux',fontsize=9)
        ax[0].grid(alpha=0.4)
        ax[0].get_yaxis().set_label_coords(-0.06,0.5)

        ax[1].plot(times_plot,targ_flux/np.mean(targ_flux),label='Target',alpha=0.7,marker='.')
        ax[1].plot(times_plot,best_ref_flux/np.mean(best_ref_flux),label='Reference',alpha=0.7,marker='.')
        ax[1].legend(bbox_to_anchor=(0.595, 0.525, 0.5, 0.5),fontsize=7)
        ax[1].set_ylabel('Norm. Flux',fontsize=9)
        ax[1].get_yaxis().set_label_coords(-0.06,0.5)

        #if use_refs == '' or use_refs == 'choice':
        for i in range(len(best_set)+1):
            if i == 0:
                lab = 'Target'
            else:
                lab = 'Ref. '+str(best_set[i-1])
            
            ax[2].plot(times_plot,x_subframe_position[i],alpha=0.5,marker='.',label=lab)
            ax[3].plot(times_plot,y_subframe_position[i],alpha=0.5,marker='.')
        ax[2].set_ylabel('X sub. pos.',fontsize=9)
        ax[2].get_yaxis().set_label_coords(-0.06,0.5)
        ax[2].legend(bbox_to_anchor=(0.5775, 0.525, 0.5, 0.5),fontsize=7)

        ax[3].set_ylabel('Y sub. pos.',fontsize=9)
        ax[3].get_yaxis().set_label_coords(-0.06,0.5)


        ax[4].plot(times_plot,targ_background,marker='.')
        ax[4].set_ylabel('Bkg.',fontsize=9)
        ax[4].get_yaxis().set_label_coords(-0.06,0.5)

        ax[5].plot(times_plot,x_seeing,color='purple',marker='.',label='X seeing')
        ax[5].plot(times_plot,y_seeing,color='orange',marker='.',alpha=0.5,label='Y seeing')
        ax[5].legend()
        ax[5].set_ylabel('Seeing (")',fontsize=9)
        ax[5].set_xlim(min(times_plot),max(times_plot))
        ax[5].get_yaxis().set_label_coords(-0.06,0.5)

        ax[6].plot(times_plot,airmass,color='gold',marker='.',alpha=0.7)
        ax[6].set_ylabel('Airmass',fontsize=9)
        ax[6].set_xlabel('Julian Date',fontsize=9)
        ax[6].set_xlim(min(times_plot),max(times_plot))
        ax[6].get_yaxis().set_label_coords(-0.06,0.5)
        plt.suptitle(target_name)
        plt.ion()
        plt.subplots_adjust(top=0.95,left=0.08,bottom=0.05,right=0.92,hspace=0.1)
        plt.savefig(local_path+'Objects/'+target_name+'/analysis/'+target_name+'_diagnostics_plot.png')
        
        # #Make plot showing zoom-ins on all nights 
        # night_locs = np.where(np.gradient(times) > 0.2)[0]
        # num_nights = int(len(night_locs)/2)+1
        # fig,ax = plt.subplots(num_nights,1,figsize=(8,10))
        # for i in range(num_nights):
        #     if i == 0:
        #         ax[i].plot(times_plot[0:night_locs[0]],best_norm_targ_flux_decorrelated[0:night_locs[0]],'k.')
        #         ax[i].set_title(target_name,fontsize=16)
        #     elif i < num_nights - 1:
        #         ax[i].plot(times_plot[night_locs[2*(i-1)+1]:night_locs[2*(i-1)+2]],best_norm_targ_flux_decorrelated[night_locs[2*(i-1)+1]:night_locs[2*(i-1)+2]],'k.')
        #     else:
        #         ax[i].plot(times_plot[night_locs[-1]:],best_norm_targ_flux_decorrelated[night_locs[-1]:],'k.')
        # plt.tight_layout()

        save_dict = {'Times':times_plot,'Flux':best_norm_targ_flux_decorrelated,'Bin starts':bin_starts,'Bin ends':bin_ends}
        pickle.dump(save_dict,open(local_path+'Objects/'+target_name+'/analysis/'+target_name+'_lightcurve.p','wb'))
        plt.show()
        pdb.set_trace()


