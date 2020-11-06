import matplotlib.pyplot as plt
plt.ion()
import pdb 
from pines_analysis_toolkit.utils import pines_dir_check, short_name_creator
from pines_analysis_toolkit.analysis.night_splitter import night_splitter
from glob import glob
from natsort import natsorted
import pandas as pd
import numpy as np

def speculoos_style_lightcurve(target, phot_type='aper'):
    '''Authors:
		Patrick Tamburo, Boston University, November 2020
	Purpose:
        Makes a lightcurve for a target using the reference star weighting technique detailed in
        Murray et al. (2020): https://arxiv.org/pdf/2005.02423.pdf
	Inputs:
        target (str): the target's full 2MASS name.
        phot_type (str, optional): 'aper' for aperture photometry
    Outputs:

	TODO:
        Make compatible with PSF photometry
    '''

    #Set up paths.
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    phot_path = pines_path/('Objects/'+short_name+'/'+phot_type+'_phot/')
    phot_files = natsorted(glob(str(phot_path/'*.csv')))
    convergence_threshold = 1e-5 

    #Loop over all photometry files.
    for i in range(len(phot_files)):
        #Read in photometry data.
        phot_file = phot_files[i]
        df = pd.read_csv(phot_file)
        df.columns = df.keys().str.strip()

        full_times = np.array(df['Time JD'])
        sources = natsorted(list(set([i.split(' ')[0]+' '+i.split(' ')[1] for i in df.keys() if (i[0] == '2') or (i[0] == 'R')])))
        num_sources = len(sources)
        num_refs = num_sources - 1

        #Split data up into individual nights. 
        night_inds = night_splitter(full_times)
        #night_inds = np.array([np.arange(len(full_times))])
        num_nights = len(night_inds)

        #Loop over each night in the dataset.
        for j in range(num_nights):
            times = full_times[night_inds[j]]
           
            #Get normalized target flux (the target's  "absolute" lightcurve).
            targ_flux = df[sources[0]+' Flux'][night_inds[j]] / np.nanmedian(df[sources[0]+' Flux'][night_inds[j]])
            
            #Read *normalized* reference fluxes into an array and generate initial weights. 
            ref_fluxes = np.zeros((num_refs, len(times)))
            ref_flux_errors = np.zeros((num_refs, len(times)))
            initial_weights = np.zeros(num_refs)
            for k in range(num_refs):
                source_name = sources[k+1]
                ref_fluxes[k] = df[source_name+' Flux'][night_inds[j]] / np.nanmedian(df[source_name+' Flux'][night_inds[j]]) #Get the reference "absolute" lightcurves.
                #initial_weights[k] = 1/(np.nanstd(ref_fluxes[k])**2) #Equation 1 of Murray et al. (2020). TODO: Modify to include other noise sources added in quadrature.
                #initial_weights[k] =1/np.mean(df[source_name+' Flux Error'][night_inds[j]])**2
                sigma =np.mean(df[source_name+' Flux Error'])/np.nanmedian(df[source_name+' Flux'][night_inds[j]])
                initial_weights[k] = 1/sigma**2

            #Normalize the initial weights such that they sum to one. 
            initial_weights = initial_weights / initial_weights.sum()
            pdb.set_trace()

            old_weights = initial_weights

            diffs = np.zeros(num_refs)+999 #Initialize to high value. 
            l = 0.
            while np.sum(diffs > convergence_threshold) > 0:
                print(l+1)
                #Create the ALC. 
                alc = np.array([np.sum(old_weights*ref_fluxes[:,i]) for i in range(len(times))]) #NOTE: should this be normalized? 

                #Generate new weights using the ALC-corrected reference lightcurves.  
                new_weights = np.zeros(num_refs)
                for k in range(num_refs):

                    corrected_flux = ref_fluxes[k] / alc
                    new_weights[k] = 1/(np.nanstd(corrected_flux)**2)

                    pdb.set_trace()


                new_weights = new_weights / new_weights.sum()
                diffs = abs(new_weights - old_weights)

                old_weights = new_weights 
                l += 1
            
            pdb.set_trace()
                