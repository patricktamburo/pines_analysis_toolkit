import matplotlib.pyplot as plt
import pdb 
from pines_analysis_toolkit.utils import pines_dir_check, short_name_creator
from glob import glob
from natsort import natsorted
import pandas as pd
import numpy as np

def raw_flux_plot(target, phot_type='aper'):
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    phot_path = pines_path/('Objects/'+short_name+'/'+phot_type+'_phot/')
    phot_files = natsorted(glob(str(phot_path/'*.csv')))

    for i in range(len(phot_files)):
        phot_file = phot_files[i]
        df = pd.read_csv(phot_file)
        df.columns = df.keys().str.strip()
        times = np.array(df['Time JD'])
        sources = natsorted(list(set([i.split(' ')[0]+' '+i.split(' ')[1] for i in df.keys() if (i[0] == '2') or (i[0] == 'R')])))
        num_sources = len(sources)

        #Unnormalized flux plot
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 5), sharex=True)
        for j in range(num_sources):
            source_name = sources[j]
            flux = df[source_name+' Flux']
            if j == 0:
                ls = '-'
            else:
                ls = ''
            ax.plot(times, flux, linestyle=ls, marker='.')
            ax.set_title('Raw flux')
        
        #Median-normalized flux plot
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 5), sharex=True)
        for j in range(num_sources):
            source_name = sources[j]
            flux = df[source_name+' Flux'] / np.median(df[source_name+' Flux'])
            if j == 0:
                ls = '-'
            else:
                ls = ''
            ax.plot(times, flux, linestyle=ls, marker='.')
            ax.set_title('Median-normalized flux')
            pdb.set_trace()

        #Mean-normalized flux plot
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 5), sharex=True)
        for j in range(num_sources):
            source_name = sources[j]
            flux = df[source_name+' Flux'] / np.mean(df[source_name+' Flux'])
            if j == 0:
                ls = '-'
            else:
                ls = ''
            ax.plot(times, flux, linestyle=ls, marker='.')
            ax.set_title('Median-normalized flux')
        pdb.set_trace()