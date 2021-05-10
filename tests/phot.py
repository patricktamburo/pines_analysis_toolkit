import pines_analysis_toolkit as pat 
import pdb
import numpy as np 


pines_path = pat.utils.pines_dir_check()
targets = ['2MASS J14392836+1929149']
#ut_dates = ['20200905', '20200906', '20200908']

for target in targets:
    short_name = pat.utils.short_name_creator(target)
    object_path = pines_path/('Objects/'+short_name)
    profile_path = object_path/(short_name.lower().replace(' ','')+'.profile')
    profile_data = pat.utils.profile_reader(profile_path)

    #Determine sources to track. 
    sources = pat.photometry.ref_star_chooser(target, profile_data['source_detect_image'], restore=True, guess_position=(profile_data['guess_position_x'], profile_data['guess_position_y']), exclude_lower_left=profile_data['exclude_lower_left'], dimness_tolerance=profile_data['dimness_tolerance'], brightness_tolerance=profile_data['brightness_tolerance'], distance_from_target=profile_data['distance_from_target'], non_linear_limit=profile_data['non_linear_limit'], edge_tolerance=profile_data['edge_tolerance'], source_detect_plot=False)

    #Find centroids of each source in every image. 
    centroided_sources = pat.photometry.centroider(target, sources, restore=False, output_plots=True, gif=False, box_w=7)

    #Generate PWV data from fyodor.
    #pat.pwv.los_pwv_data_generator(target, ut_dates, centroided_sources)

    #Perform photometry using a variety of methods. 
    pat.photometry.aper_phot.fixed_aper_phot(target, centroided_sources, ap_radii=[4.0, 4.5, 5.0, 5.5, 6.0], an_in=9, an_out=30)
    #pat.photometry.aper_phot.variable_aper_phot(target, centroided_sources, [1.0])
    #pat.photometry.epsf_phot(target, centroided_sources, plots=True)
    #pat.photometry.basic_psf_phot(target, centroided_sources, plots=True)

    # #Create a lightcurve. 
    pat.analysis.weighted_lightcurve(target, phot_type='aper', mode='night', plots=True)
    pat.analysis.weighted_lightcurve(target, phot_type='aper', mode='global', plots=True)

    # #Generate diagnostic plots. 
    pat.analysis.diagnostic_plots.relative_cutout_position_plot(target, centroided_sources)
    pat.analysis.diagnostic_plots.absolute_image_position_plot(target, centroided_sources)
    pat.analysis.diagnostic_plots.background_plot(target, centroided_sources) #TODO FIX reading photometry in. 
    pat.analysis.diagnostic_plots.seeing_plot(target, centroided_sources)
    pat.analysis.diagnostic_plots.airmass_plot(target, centroided_sources)
    
    #Create plots of corrected reference star fluxes. 
    pat.analysis.analysis_plots.corr_all_sources_plot(target)
    
    #Generate DV report. 
    pat.output.dv_report(target)