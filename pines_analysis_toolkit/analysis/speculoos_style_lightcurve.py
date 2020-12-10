import matplotlib.pyplot as plt
colormap = plt.cm.viridis
import pdb 
from pines_analysis_toolkit.utils import pines_dir_check, short_name_creator
from pines_analysis_toolkit.analysis.night_splitter import night_splitter
from glob import glob
from natsort import natsorted
import pandas as pd
import numpy as np
import time 
import julian 
import matplotlib.dates as mdates
import os
import shutil

def speculoos_style_lightcurve(target, sources, output_plots=False, phot_type='aper'):
    '''Authors:
		Patrick Tamburo, Boston University, November 2020
	Purpose:
        Makes a lightcurve for a target using the reference star weighting technique detailed in
        Murray et al. (2020): https://arxiv.org/pdf/2005.02423.pdf
	Inputs:
        target (str): the target's full 2MASS name.
        sources (pandas dataframe): df of source names and positions, for weighting based on distances. 
        phot_type (str, optional): 'aper' for aperture photometry.
    Outputs:

	TODO:
        Make compatible with PSF photometry
    '''

    def absolute_plot():
        '''Makes a plot of comparison star absolute lightcurves along with ALC.'''
        date = df['Time UT'][night_inds[j]][0].split('T')[0].strip()
        fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(15,6), sharex=True)
        for k in range(len(absolute_ref_fluxes)):
            lab = '{:6s} - {:4.3f}'.format('CS '+str(k+1), old_weights[k])
            color = colormap((k+1)*int(256/len(absolute_ref_fluxes-1))-1)
            ax.plot(times, absolute_ref_fluxes[k], linestyle='', marker='.', label=lab, color=color, ms=((old_weights[k]/np.max(old_weights))**0.25)*15)
        ax.plot(times, alc, lw=3, marker='s', label='ALC', color='k', ms=8, markeredgewidth=1.5, markerfacecolor='none', linestyle='')
        ax.set_ylabel('Absolute Flux', fontsize=14)
        ax.set_xlabel('Time', fontsize=14)
        ax.legend(fancybox=True)
        ax.set_title(source_names[0]+', '+date+'\nAbsolute Lightcurves and ALC, Iteration '+str(int(l)), fontsize=16)
        myFmt = mdates.DateFormatter('%H:%M')
        ax.xaxis.set_major_formatter(myFmt)
        fig.autofmt_xdate()
        ax.tick_params(labelsize=12)
        plt.tight_layout()
        filename = pines_path/('Objects/'+short_name+'/analysis/speculoos_absolute_plots/'+str(int(l)).zfill(3)+'.png')
        plt.savefig(filename)
        plt.close()

    def differential_plot():
        '''Makes a plot of each comparison star's differential lightcurve (absolute corrected by ALC).'''
        date = df['Time UT'][night_inds[j]][0].split('T')[0].strip()
        fig,ax = plt.subplots(nrows=num_sources, ncols=1, figsize=(15,20), sharex=True, sharey=True)
        ax[0].plot(times, targ_flux/alc, marker='.', ms=10, label='Target', color='r', linestyle='')
        ax[0].axhline(y=1, color='k', zorder=0)
        ax[0].legend(fancybox=True, loc='upper right', fontsize=14)
        for k in range(num_refs):
            color = colormap((k+1)*int(256/len(absolute_ref_fluxes-1))-1)
            lab = 'CS '+str(k+1)
            ax[k+1].plot(times, differential_ref_fluxes[k], linestyle='', marker='.', color=color, label=lab, ms=10)
            ax[k+1].tick_params(labelsize=12)
            ax[k+1].legend(fancybox=True, loc='upper right', fontsize=14)
            ax[k+1].axhline(y=1, color='k', zorder=0)
            ax[k+1].text(1.02, 0.5, "$\sigma$ = {:1.3f}".format(np.std(differential_ref_fluxes[k])), transform=ax[k+1].transAxes, fontsize=14, bbox=dict(facecolor=color, alpha=0.5)).set_clip_on(False)
        ax[int((num_refs-1)/2)].set_ylabel('Differential Flux', fontsize=14)
        ax[num_refs-1].set_xlabel('Time (UT)', fontsize=14)
        ax[0].set_title(source_names[0]+', '+date+'\nDifferential Lightcurves, Iteration '+str(int(l)), fontsize=16)
        myFmt = mdates.DateFormatter('%H:%M')
        ax[num_refs-1].xaxis.set_major_formatter(myFmt)
        fig.autofmt_xdate()
        plt.tight_layout()
        fig.subplots_adjust(hspace=0.1, right=0.9)
        filename = pines_path/('Objects/'+short_name+'/analysis/speculoos_differential_plots/'+str(int(l)).zfill(3)+'.png')
        plt.savefig(filename)
        plt.close()
    
    def pretty_printer(init=False):
        '''Prints out comparison star weights/differences in table format. '''
        if init:
            print_str = ''
            print_str_2 = ''
            for k in range(num_refs):
                if k == 0:
                    print_str += '| Source |  CS'+str(k+1)+'   |'
                else:
                    print_str += '  CS'+str(k+1)+'   |'
                if k+1 >= 10:
                    print_str_2 += '  {:>1.2f}   |'.format(initial_weights[k])
                else:
                    if k == 0:
                        print_str_2 += '| Weight |  {:>1.2f}  |'.format(initial_weights[k])

                    else:
                        print_str_2 += '  {:>1.2f}  |'.format(initial_weights[k])
            
            print('INITIAL WEIGHTS ')
            print('-'*len(print_str))
            print(print_str)
            print('-'*len(print_str))
            print(print_str_2)
            print('-'*len(print_str))
            print('')
            time.sleep(1)

        else:
            print_str = ''
            print_str_2 = ''
            print_str_3 = ''
            for k in range(num_refs):
                if k == 0:
                    print_str += '| Source |  CS'+str(k+1)+'   |'
                else:
                    print_str += '  CS'+str(k+1)+'   |'

                if k+1 >= 10:
                    print_str_2 += '  {:>1.2f}   |'.format(new_weights[k])
                    diff = new_weights[k]-old_weights[k]
                    if diff >= 0:
                        print_str_3 += ' +{:>1.2f}   |'.format(diff)
                    else:
                        print_str_3 += ' {:>1.2f}   |'.format(diff)
                else:
                    if k == 0:
                        print_str_2 += '| Weight |  {:>1.2f}  |'.format(new_weights[k])
                        diff = new_weights[k]-old_weights[k]
                        if diff >= 0:
                            print_str_3 += '| Diff   | +{:>1.2f}  |'.format(diff)
                        else:
                            print_str_3 += '| Diff   | {:>1.2f}  |'.format(diff)
                    else:
                        print_str_2 += '  {:>1.2f}  |'.format(new_weights[k])
                        diff = new_weights[k]-old_weights[k]
                        if diff >= 0:
                            print_str_3 += ' +{:>1.2f}  |'.format(diff)
                        else:
                            print_str_3 += ' {:>1.2f}  |'.format(diff)

            print('ITERATION '+str(int(l+1)))
            print('-'*len(print_str))
            print(print_str)
            print('-'*len(print_str))
            print(print_str_2)
            print(print_str_3)
            print('-'*len(print_str))
            print('')

    #Set up paths and get photometry files. 
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    phot_path = pines_path/('Objects/'+short_name+'/'+phot_type+'_phot/')
    phot_files = natsorted(glob(str(phot_path/'*.csv')))

    #Set up directories for output plots. 
    if output_plots:
        subdirs = glob(str(pines_path/('Objects/'+short_name+'/analysis'))+'/speculoos*/')
        #Delete any source directories that are already there. 
        for name in subdirs:
            shutil.rmtree(name)
        #Create new source directories.
        os.mkdir(pines_path/('Objects/'+short_name+'/analysis/speculoos_absolute_plots/'))
        os.mkdir(pines_path/('Objects/'+short_name+'/analysis/speculoos_differential_plots/'))
        
    convergence_threshold = 1e-5 #Following Murray et al. (2020). 

    #Get pixel distances between comparison stars and target for weighting purposes. 
    targ_x = sources['Source Detect X'][0]
    targ_y = sources['Source Detect Y'][0]
    pixel_dists = np.array(np.sqrt((targ_x - sources['Source Detect X'][1:])**2 + (targ_y - sources['Source Detect Y'][1:])**2))

    #Loop over all photometry files.
    for i in range(len(phot_files)):
        #Read in photometry data.
        phot_file = phot_files[i]
        df = pd.read_csv(phot_file)
        df.columns = df.keys().str.strip()
        full_times = np.array(df['Time JD'])
        source_names = natsorted(list(set([i.split(' ')[0]+' '+i.split(' ')[1] for i in df.keys() if (i[0] == '2') or (i[0] == 'R')])))
        ref_names = source_names[1:]
        num_sources = len(source_names)
        num_refs = len(ref_names)

        #Split data up into individual nights. 
        night_inds = night_splitter(full_times)
        num_nights = len(night_inds)

        #Loop over each night in the dataset.
        for j in range(num_nights):
            times = full_times[night_inds[j]]
            times = np.array([julian.from_jd(times[i], fmt='jd') for i in range(len(times))])

            #Get normalized target flux (the target's  "absolute" lightcurve).                
            t_flux = np.array([float(i) for i in df[source_names[0]+' Flux'][night_inds[j]]])
            targ_flux = t_flux / np.nanmean(t_flux)
            
            #Read *normalized* reference fluxes into an array and generate initial variability weights (w_var). 
            raw_ref_fluxes = np.zeros((num_refs, len(times)))
            absolute_ref_fluxes = np.zeros((num_refs, len(times)))
            initial_w_var = np.zeros(num_refs)
            for k in range(num_refs):
                source_name = ref_names[k]
                raw_ref_fluxes[k] = np.array([float(i) for i in df[source_name+' Flux'][night_inds[j]]]) #Get the raw reference lightcurve.
                norm = np.nanmean(raw_ref_fluxes[k])
                absolute_ref_fluxes[k] = raw_ref_fluxes[k] / norm #Get the reference "absolute" lightcurve.
                err = [float(i) for i in df[source_name+' Flux Error'][night_inds[j]]] #Get the error on the reference lightcurve
                sigma = np.nanmean(err)/norm
                initial_w_var[k] = 1/sigma**2 #Set initial variability weights using the error on the reference lightcurve. 
                plt.errorbar(times, raw_ref_fluxes[k], err, marker='.', linestyle='')
                plt.yscale('log')

            pdb.set_trace()
            #Create initial weights based on comparison star distances from the target (w_dist).
            initial_w_dist = 1 / (1 + (pixel_dists/np.max(pixel_dists))**2)

            initial_w_dist = initial_w_dist / initial_w_dist.sum()
            
            #Combine w_var and w_dist into one weight. 
            initial_weights = initial_w_var * initial_w_dist

            #Normalize the initial weights such that they sum to one. 
            initial_weights = initial_weights / initial_weights.sum()
            pretty_printer(init=True)

            #Perform an iterative loop until weights converge:
            #   1. Create the artificial comparison lightcurve (ALC) from the weighted mean of the absolute (normalized) flux of the n comparison stars in each frame. 
            #   2. Divide every comparison star's absolute lightcurve by this ALC to produce a differential lightcurve. 
            #   3. Replace w_var for each comparison star using the *measured* standard deviation in the differential lightcurves: w'_var,i = 1/sigma_i^2.
            #   4. Repeat until weights are constant to within threshold value of 1e-5. 

            #Set up stuff for the while loop.s
            old_weights = initial_weights
            diffs = np.zeros(num_refs)+999 #Initialize to high value. 
            l = 0.
            while np.sum(diffs > convergence_threshold) > 0:
                #Create the ALC. 
                alc = np.array([np.sum(old_weights*absolute_ref_fluxes[:,i]) for i in range(len(times))]) 

                #Generate new weights using the ALC-corrected reference lightcurves.  
                new_weights = np.zeros(num_refs)
                sigmas = np.zeros(num_refs)
                differential_ref_fluxes = np.zeros((num_refs, len(times)))

                for k in range(num_refs):
                    differential_ref_fluxes[k] = absolute_ref_fluxes[k] / alc #Create differential comparison star lightcurves.         
                    sigmas[k] = np.std(differential_ref_fluxes[k]) #Measure the standard deviation of the differential lightcurves. 
                    new_weights[k] = 1/(sigmas[k]**2) * initial_w_dist[k] #Update the weights using the measured standard deviations and the original distance weights. 
               
                if output_plots: #Can make output plots of the absolute lightcurves + ALC, differential lightcurves. 
                    absolute_plot()
                    differential_plot()
                
                new_weights = new_weights / new_weights.sum() #Normalize the new weights to sum to 1. 

                pretty_printer()

                diffs = abs(new_weights - old_weights) #Calculate the difference between the new weights and the old weights. Want max(diffs) < threshold for loop to end. 
                old_weights = new_weights #Update the old weights and start the loop over. 
                
                l += 1 

            weights = new_weights
            print('ALC weights converged to within {}.'.format(convergence_threshold))
            pdb.set_trace()

    