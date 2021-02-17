from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
import pdb 
import shutil 

def output_wrangler(target):
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    object_path = pines_path/('Objects/'+short_name)
    sources_path = object_path/('sources')
    aper_phot_path = object_path/('aper_phot')
    analysis_path = object_path/('analysis')
    output_path = object_path/('output')

    #Copy souce detect image. 
    src = sources_path/('target_and_refs.png')
    dest = output_path/('target_and_refs.png')
    shutil.copyfile(src, dest)

    #Copy image centroids. 
    src = sources_path/('target_and_references_centroids.csv')
    dest = output_path/('target_and_references_centroids.csv')
    shutil.copyfile(src, dest)

    #Copy best aperture raw photometry.
    best_ap_path = analysis_path/('optimal_aperture.txt')
    with open(best_ap_path, 'r') as f:
        best_ap = f.readlines()[0]
    src = aper_phot_path/(short_name+'_aper_phot_'+best_ap+'_pix_radius.csv')
    dest = output_path/(short_name+'_aper_phot_'+best_ap+'_pix_radius.csv')
    shutil.copyfile(src, dest)

    #Copy best aperture weighted lightcurve photometry. 
    src = analysis_path/('aper_phot_analysis/'+best_ap+'/'+short_name+'_weighted_lc_aper_phot_'+best_ap+'_pix.csv')
    dest = output_path/(short_name+'_weighted_lc_aper_phot_'+best_ap+'_pix.csv')
    shutil.copyfile(src, dest)

    pdb.set_trace()
    return