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
    #pat.photometry.aper_phot.fixed_aper_phot(target, centroided_sources, [3.5, 4.5], an_in=9)
    pat.photometry.aper_phot.variable_aper_phot(target, centroided_sources)
    #pat.photometry.epsf_phot(target, centroided_sources, plots=True)
    #pat.photometry.basic_psf_phot(target, centroided_sources, plots=True)

    #Create a lightcurve. 
    #pat.analysis.simple_lightcurve(target, sources, centroided_sources, plot_mode='separate')
    #pat.analysis.weighted_lightcurve(target, phot_type='aper', mode='global', plots=True)

    #Generate diagnostic plots. 
    # pat.analysis.diagnostic_plots.seeing_plot(target, centroided_sources)
    # pat.analysis.diagnostic_plots.relative_cutout_position_plot(target, centroided_sources)
    # pat.analysis.diagnostic_plots.absolute_image_position_plot(target, centroided_sources)
    # pat.analysis.diagnostic_plots.background_plot(target, centroided_sources)

    #Generate DV report. 
    #pat.output.dv_report(target)

    #Move all final data products to the target's output directory. 
    pat.output.output_wrangler(target)
