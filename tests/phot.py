from pines_analysis_toolkit import output
import pines_analysis_toolkit as pat 
from pathlib import Path
import pdb 

target = '2MASS J16202614-0416315'

pines_path = pat.utils.pines_dir_check()
short_name = pat.utils.short_name_creator(target)
profile_path = pines_path/('Objects/'+short_name+'/'+short_name.replace(' ','').lower()+'.profile')
profile_data = pat.utils.profile_reader(profile_path)

#Get reference stars
sources = pat.photometry.ref_star_chooser(target, profile_data['source_detect_image'], restore=True, guess_position=(profile_data['guess_position_x'],profile_data['guess_position_y']), non_linear_limit=profile_data['non_linear_limit'], dimness_tolerance=profile_data['dimness_tolerance'], brightness_tolerance=profile_data['brightness_tolerance'], distance_from_target=profile_data['distance_from_target'], edge_tolerance=profile_data['edge_tolerance'], exclude_lower_left=profile_data['exclude_lower_left'])

#Get an astrometric solution for the source detection image so that we can figure out the spectral types of our reference stars. 
api_key = 'vqrhxpyypnpmqbzo'
source_detect_image_path = Path(pines_path/('Objects/'+short_name+'/reduced/'+profile_data['source_detect_image']))
pat.astrometry.source_detect_astrometry(target, api_key, source_detect_image_path)

#Get world coordinates of sources in the source_detect_image. 
sources = pat.astrometry.source_pixels_to_world(target, source_detect_image_path)

#Query 2MASS and Gaia for information on sources. 
sources = pat.photometry.gaia_cmd(target, sources, red_only=False)

#Get centroids for sources.
centroided_sources = pat.photometry.centroider(target, sources, restore=True, output_plots=False)

#Do photometry on centroided sources.
pat.photometry.aper_phot.fixed_aper_phot(target, centroided_sources, ap_radii=[4.0, 4.5, 5.0, 5.5, 6.0])

#Make lightcurves. 
pat.analysis.weighted_lightcurve(target, plots=True)
pdb.set_trace()