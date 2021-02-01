import pines_analysis_toolkit as pat
from pdb import set_trace as stop

targets = ['2MASS J06420559+4101599']
target = targets[0]
dates = ['20201005', '20201007'] #UT dates the targets were observed on.


exptimes = [60.] #Exposure time of the images in seconds. 
flat_date = '20201003' #The date that dome flat calibration data were taken for the target. 
dark_date = '20201003' #The date that dark calibration data were taken for the target.
band = 'J'

#sftp = pat.utils.pines_login()

# #for exptime in exptimes:
# #    pat.data.dark(dark_date, exptime, upload=False, delete_raw=True, sftp=sftp, dark_start=1, dark_stop=10)
# pat.data.dome_flat_field(flat_date, band, sftp=sftp, upload=True, delete_raw=True, lights_on_start=990, lights_on_stop=1089, lights_off_start=1090, lights_off_stop=1189)

# for exptime in exptimes:
#    pat.data.variable_pixels(dark_date, exptime, upload=True, sftp=sftp)
#    pat.data.hot_pixels(dark_date, exptime, box_l=5, upload=True, sftp=sftp)

# pat.data.dead_pixels(flat_date, band, upload=False, sftp=sftp)

# for exptime in exptimes:
#    pat.data.bpm_maker(flat_date, dark_date, exptime, band, upload=True, sftp=sftp)

# for target in targets:
#     pat.data.get_raw_science_files(sftp, target)
#     pat.data.reduce(target, delete_raw=True, delete_reduced=False, upload=True, sftp=sftp)

# sftp.close()


#target = targets[0]

# #Update logs with more accurate shifts.
# for target in targets:
#     for date in dates:
#         pat.observing.log_updater(target, date, upload=True)

# stop()

sources = pat.photometry.ref_star_chooser(target, restore=True, source_detect_image_ind=30, 
            exclude_lower_left=False, dimness_tolerance=0.9, distance_from_target=300., non_linear_limit=3300, 
            edge_tolerance=80., source_detect_plot=False)


centroided_sources = pat.photometry.centroider(target, sources, restore=True, output_plots=False, gif=False, box_w=7)

#pat.photometry.aper_phot(target, centroided_sources, [4], an_in=9)
#pat.photometry.epsf_phot(target, centroided_sources, plots=True)
#pat.photometry.basic_psf_phot(target, centroided_sources, plots=True)

#pat.analysis.lightcurve(target, sources, centroided_sources)
#pat.analysis.speculoos_style_lightcurve(target, sources, output_plots=True, phot_type='aper')

for target in targets:
    print(target)
    pat.analysis.simple_lightcurve(target, sources, centroided_sources, plot_mode='separate')
    pat.analysis.weighted_lightcurve(target, phot_type='aper', mode='night')
