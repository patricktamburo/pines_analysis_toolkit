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

def noise_estimator():

    target_name = '2MASS 0344+0111'
    G = 8.21 #e- / DN
    #RN = 19 #e- / pix. From Clemens (2007).
    RN = 32.89600862127023 #e- / pix. Measured on a set of biases from Nov. 2019.
    t = 30 #s
    r = 3 #pix
    npix = np.pi*r**2 #pix
    D = 0.98 #e- / pix / s

    lightcurve_std = 0.02265514985623077
    #Pathing
    path = 'P:\\Photometry\\PINES\\'+target_name+'\\aper_phot\\photometry\\'+target_name+'_ap='+str(r)+'_photometry.p'
    analysis_path = 'P:\\Photometry\\PINES\\'+target_name+'\\analysis\\'
    restored_output = pickle.load(open(path,'rb'))
    targ_flux = restored_output['Flux'][0]
    background = restored_output['Background'][0]
    times = restored_output['Time']
    airmass = restored_output['Airmass']
    avg,med,std = sigma_clipped_stats(targ_flux)
    #pdb.set_trace()

    use_locs_2 = np.where( (airmass < 2) & (background < 6000) & (targ_flux > (med-3*std)))[0] #For 2mass 1941 late sep 2019

    targ_flux = targ_flux[use_locs_2]
    background = background[use_locs_2]
    times = times[use_locs_2]

    #Read in lightcurve_data
    lc = pickle.load(open(analysis_path+target_name+'_lightcurve.p','rb'))
    lc_times = lc['Times']
    lc_flux = lc['Flux']
    lc_bin_starts = np.array([int(i) for i in lc['Bin starts']])
    lc_bin_ends = np.array([int(i) for i in lc['Bin ends']])
    #Get the date of each image
    days = [lc_times[i].day for i in range(len(lc_times))]
    unique_days = np.unique(days)

    R_star = targ_flux * G / t #Convert to e-/s
    R_sky = background * G / t #Convert to e-/s
    pdb.set_trace()

    source_photon_noise = np.sqrt(R_star*t)
    sky_photon_noise = np.sqrt(npix*R_sky*t)
    read_noise = np.zeros(len(source_photon_noise))+np.sqrt(npix*RN**2)
    dark_current = np.zeros(len(source_photon_noise))+np.sqrt(D*npix*t)

    SNR = R_star*t/np.sqrt((R_star*t)+(R_sky*t*npix)+(npix*RN**2)+(D*npix*t))
    noise = 1/SNR
    star_component = 1/(R_star*t)
    sky_component = R_sky*npix*t/(R_star*t)**2
    rn_component = npix*RN**2/(R_star*t)**2
    dark_component = D*npix*t/(R_star*t)**2

    plt.ion()
    for i in range(len(unique_days)):
        date = unique_days[i]
        where_date = np.where(days == date)[0] 
        title = str(lc_times[where_date[0]].month)+'-'+str(lc_times[where_date[0]].day)+'-'+str(lc_times[where_date[0]].year)
        lc_variance = np.std(lc_flux[where_date])**2
        
        fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(9,6))
        ax.plot(times[where_date],noise[where_date]**2,lw=2,label='Calculated variance',zorder=1,marker='o')
        ax.plot(times[where_date],sky_component[where_date],'k',label='Sky',zorder=0,marker='o')
        ax.plot(times[where_date],star_component[where_date],'b',label='Star',zorder=0,marker='o')
        ax.plot(times[where_date],rn_component[where_date],'r',label='Read noise',zorder=0,marker='o')
        ax.plot(times[where_date],dark_component[where_date],'m',label='Dark current',zorder=0,marker='o')
        ax.plot([times[where_date[0]],times[where_date[-1]]],[lc_variance,lc_variance],color='orange',lw=2,zorder=2,label='Measured variance')
        ax.set_title(title)
        plt.legend()
        plt.ylabel('$\sigma^2$',fontsize=14)
        plt.xlabel('Time',fontsize=14)
        plt.tight_layout()


    plt.show()
