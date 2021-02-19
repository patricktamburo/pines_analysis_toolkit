import pines_analysis_toolkit as pat
import pdb 

targets = ['2MASS J08433323+1024470']
dates = ['20210204'] #UT dates the targets were observed on.
exptimes = [30., 60.] #Exposure time of the images in seconds. 
flat_date = '20210204' #The date that dome flat calibration data were taken for the target. 
dark_date = '20210204' #The date that dark calibration data were taken for the target.
band = 'J'

#sftp = pat.utils.pines_login()

# for exptime in exptimes:
#     pat.data.dark(dark_date, exptime, upload=True, delete_raw=True, sftp=sftp)
# pat.data.dome_flat_field(flat_date, band, sftp=sftp, upload=True, delete_raw=True)

# for exptime in exptimes:
#    pat.data.variable_pixels(dark_date, exptime, upload=True, sftp=sftp)
#    pat.data.hot_pixels(dark_date, exptime, box_l=5, upload=True, sftp=sftp)

# pat.data.dead_pixels(flat_date, band, upload=False, sftp=sftp)

# for exptime in exptimes:
#    pat.data.bpm_maker(flat_date, dark_date, exptime, band, upload=True, sftp=sftp)

# for target in targets:
#     pat.data.get_raw_science_files(sftp, target)
#     pat.data.reduce(target, delete_raw=True, delete_reduced=False, upload=False, sftp=sftp)

# pdb.set_trace()

#Update logs with more accurate shifts.
#for target in targets:
#   for date in dates:
#       pat.observing.log_updater(target, date, sftp, upload=True)

#sftp.close()

# pdb.set_trace()

for target in targets:

    #Determine sources to track. 
    sources = pat.photometry.ref_star_chooser(target, restore=True, source_detect_image_ind=30, 
                exclude_lower_left=False, dimness_tolerance=0.75, distance_from_target=500., non_linear_limit=3300, 
                edge_tolerance=40., source_detect_plot=False)

    #Find centroids of each source in every image. 
    centroided_sources = pat.photometry.centroider(target, sources, restore=True, output_plots=False, gif=False, box_w=7)

    #Perform photometry using a variety of methods. 
    #pat.photometry.aper_phot.fixed_aper_phot(target, centroided_sources, [4.2, 4.3, 4.4], an_in=9, an_out=30)
    #pat.photometry.aper_phot.variable_aper_phot(target, centroided_sources, [0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3])
    #pat.photometry.epsf_phot(target, centroided_sources, plots=True)
    #pat.photometry.basic_psf_phot(target, centroided_sources, plots=True)

    #Create a lightcurve. 
    #pat.analysis.simple_lightcurve(target, sources, centroided_sources, plot_mode='separate')
    #pat.analysis.weighted_lightcurve(target, phot_type='aper', mode='night', plots=False)
    #pat.analysis.weighted_lightcurve(target, phot_type='aper', mode='global', plots=False)

    #Generate diagnostic plots. 
    #pat.analysis.diagnostic_plots.relative_cutout_position_plot(target, centroided_sources)
    #pat.analysis.diagnostic_plots.absolute_image_position_plot(target, centroided_sources)
    #pat.analysis.diagnostic_plots.background_plot(target, centroided_sources) #TODO FIX reading photometry in. 
    #pat.analysis.diagnostic_plots.seeing_plot(target, centroided_sources)
    #pat.analysis.diagnostic_plots.airmass_plot(target, centroided_sources)
    
    #Create plots of corrected reference star fluxes. 
    #pat.analysis.analysis_plots.corr_all_sources_plot(target)
    
    #Generate DV report. 
    pat.output.dv_report(target)

    #Move all final data products to the target's output directory. 
    #pat.output.output_wrangler(target)
