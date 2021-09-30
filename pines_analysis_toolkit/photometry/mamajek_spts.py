import pines_analysis_toolkit as pat
import pandas as pd 
import numpy as np 
def mamajek_spts(sources):
    pines_path = pat.utils.pines_dir_check()
    mamajek_path = pines_path/('Misc/mamajek_spts.txt')
    mamajek_df = pd.read_csv(mamajek_path, sep=r"[ ]{1,}")
    mamajek_spectral_types = np.array(mamajek_df['#SpT'])
    mamajek_teffs = np.array(mamajek_df['Teff'])

    mamajek_bp_rp = np.array(mamajek_df['Bp-Rp'])
    mamajek_bp_rp[mamajek_bp_rp  == '...'] = np.nan
    mamajek_bp_rp = np.array([float(i) for i in mamajek_bp_rp])

    source_bp_rps = np.array(sources['bp_rp'])
    source_spts = []
    source_teffs = []
    for i in range(len(source_bp_rps)):
        source_bp_rp = source_bp_rps[i] 
        if (not np.isnan(source_bp_rp)) and (source_bp_rp > -0.120) and (source_bp_rp < 4.0):
            closest_spt_ind = np.where(abs(mamajek_bp_rp-source_bp_rp) == np.nanmin(abs(mamajek_bp_rp-source_bp_rp)))[0][0]
            source_spts.append(mamajek_spectral_types[closest_spt_ind])
            source_teffs.append(mamajek_teffs[closest_spt_ind])
        else:
            source_spts.append(np.nan)
            source_teffs.append(np.nan)
    
    sources['SpT'] = source_spts
    sources['Teff'] = source_teffs
    
    return sources