import pdb
from pines_analysis_toolkit.photometry.mamajek_spts import mamajek_spts 
import pines_analysis_toolkit as pat 
from astroquery.vizier import Vizier
from astroquery.gaia import Gaia
import astropy.units as u 
from astropy.coordinates import SkyCoord, Angle
import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd
plt.ioff()

def gaia_cmd(target, sources, catalog='eDR3', plot=True, force_output_path=''):
    '''Authors:
            Patrick Tamburo, Boston University, May 2021
        Purpose:
            Queries Gaia for sources and creates a CMD for them. 
            Adds Gaia M_G and BP-RP to the source csv file. 
        Inputs:
            target (str): The target's full 2MASS name.
            sources (pd.DataFrame): DataFrame of sources to be tracked in the field. These sources must have RA and Dec values.
            catalog (str, optional): The name of the Gaia catalog you wish to query, 'DR2' or 'eDR3'.
            plot (bool, optional): Whether or not to save the CMD to the object's 'sources' directory.
            red_only (bool, optional): Whether or not to limit the retained sources to the reddest stars in the field. 
            force_output_path (path): if you want to manually set an output directory, specify the top-level here (i.e. the directory containing Logs, Objects, etc.)
        Outputs:
            sources (pd.DataFrame): DataFrame of sources, with M_G and BP-RP columns added. Also updates the target_and references_source_detection.csv file in the object's 'sources' directory. 
        TODO:
            None.
    '''

    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pat.utils.pines_dir_check()

    short_name = pat.utils.short_name_creator(target)

    #Set the Gaia data release to query. 
    if catalog == 'DR2':
        Gaia.MAIN_GAIA_TABLE = "gaiadr2.gaia_source"  
    elif catalog == 'eDR3':
        Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_source" 
    else:
        raise ValueError('catalog must be DR2 or eDR3.')

    num_s = len(sources)
    bp_rp = np.zeros(num_s) #Gaia Bp - Rp color
    p_mas = np.zeros(num_s) #Parallax in mas
    M_G = np.zeros(num_s) #Absolute Gaia Rp 

    print('Querying sources in Gaia {}.'.format(catalog))

    for i in range(num_s):
        ra = sources['Source Detect RA'][i]
        dec = sources['Source Detect Dec'][i]
        c = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))

        #Query Gaia 
        #query_gaia = Gaia.cone_search(c, radius=5*u.arcsec)
        #result_gaia = query_gaia.get_results()
        
        result_gaia = Gaia.query_object_async(coordinate=c, width=u.Quantity(5, u.arcsec), height=u.Quantity(5, u.arcsec))

        #If no sources at this position, return NaNs.
        if len(result_gaia) == 0:
            bp_rp[i] = np.nan
            p_mas[i] = np.nan
            M_G[i] = np.nan
            print('No source found in Gaia {} at ({:1.4f},{:1.4f}). Returning NaNs.'.format(catalog, ra, dec))
        #If one source found, return values
        else:
            bp_rp[i] = result_gaia['bp_rp'][0]
            m_G = result_gaia['phot_g_mean_mag'][0]
            p_mas[i] = result_gaia['parallax'][0]
            dist_pc = 1/(p_mas[i]/1000)
            M_G[i] = m_G + 5 - 5*np.log10(dist_pc) 

    if plot:
        plt.figure(figsize=(6,9))
        plt.scatter(bp_rp, M_G)
        plt.title(sources['Name'][0]+' field CMD', fontsize=18)
        plt.xlabel('G$_{BP}$ - G$_{RP}$', fontsize=16)
        plt.ylabel('M$_G$', fontsize=16)
        plt.xlim(-1,5)
        plt.ylim(-5,16)
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(pines_path/('Objects/'+short_name+'/sources/gaia_cmd.png'), dpi=300)

    sources['M_G'] = M_G
    sources['bp_rp'] = bp_rp

    sources = mamajek_spts(sources)
    
    sources.to_csv(pines_path/('Objects/'+short_name+'/sources/target_and_references_source_detection.csv'), index=0, na_rep='NaN')

    return sources