import pines_analysis_toolkit as pat
import pdb

#targets = ['2MASS J01365662+0933473']
targets = ['SIMP 0136'] #Declare a target. This is SIMP J013656.5+093347.
exptimes = [30.] #Exposure time of the images in seconds. 
flat_date = '20151107' #The date that dome flat calibration data were taken for the target. 
dark_date = '20151110' #The date that dark calibration data were taken for the target.
band = 'J'

#sftp = pat.utils.pines_login()

#Do calibration stuff. 
#for exptime in exptimes:
    #  pat.data.dark(dark_date, exptime, dark_start=1966, dark_stop=1994, upload=False, delete_raw=False, sftp=sftp)
# pat.data.dome_flat_field(flat_date, band, lights_on_start=4836, lights_on_stop=4935, lights_off_start=4936, lights_off_stop=5035, sftp=sftp, upload=False, delete_raw=True)

# for exptime in exptimes:
#     pat.data.variable_pixels(dark_date, exptime, upload=False, sftp=sftp)
#     pat.data.hot_pixels(dark_date, exptime, box_l=5, upload=False, sftp=sftp)

# pat.data.dead_pixels(flat_date, band, upload=False, sftp=sftp)

# for exptime in exptimes:
#     pat.data.bpm_maker(flat_date, dark_date, exptime, band, upload=False, sftp=sftp)

# for target in targets:
#  pat.data.get_raw_science_files(sftp, target)
#  pat.data.reduce(target, delete_raw=True, delete_reduced=False, upload=False)

# sftp.close()

# pdb.set_trace()

target = targets[0]
sources = pat.photometry.ref_star_chooser(target, restore=False, source_detect_image_ind=30, exclude_lower_left=False, dimness_tolerance=0.25, distance_from_target=900., non_linear_limit=3500, edge_tolerance=100.)

centroided_sources = pat.photometry.centroider(target, sources, restore=False, output_plots=True, gif=False, box_w=7)

pat.photometry.aper_phot(target, centroided_sources, [5], an_in=9)
#pat.photometry.epsf_phot(target, centroided_sources, plots=True)
#pat.photometry.basic_psf_phot(target, centroided_sources, plots=True)

#pat.analysis.lightcurve(target, sources, centroided_sources)
#pat.analysis.speculoos_style_lightcurve(target, sources, output_plots=True, phot_type='aper')
pat.analysis.pines_alc(target, phot_type='aper', mode='night')