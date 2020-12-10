import pines_analysis_toolkit as pat
import pdb

#targets = ['2MASS J01365662+0933473']
targets = ['2MASS J05012406-0010452']
exptimes = [60.]
band = 'J'
flat_date = '20201203'
dark_date = '20201203'

#sftp = pat.utils.pines_login()

#Do calibration stuff. 
# for exptime in exptimes:
#   pat.data.dark(dark_date, exptime, upload=True, delete_raw=True, sftp=sftp)
#pat.data.dome_flat_field(flat_date, band, sftp=sftp, upload=True, delete_raw=True)

# for exptime in exptimes:
#     pat.data.variable_pixels(dark_date, exptime, upload=True, sftp=sftp)
#     pat.data.hot_pixels(dark_date, exptime, box_l=5, upload=True, sftp=sftp)

#pat.data.dead_pixels(flat_date, band, upload=True, sftp=sftp)

# for exptime in exptimes:
#   pat.data.bpm_maker(flat_date, dark_date, exptime, band, upload=True, sftp=sftp)

# for target in targets:
#  pat.data.get_raw_science_files(sftp, target)
#  pat.data.reduce(target, delete_raw=False, delete_reduced=False, upload=False, sftp=sftp)

# sftp.close()

#pdb.set_trace()

target = targets[0]
sources = pat.photometry.ref_star_chooser(target, restore=False, source_detect_image_ind=60, exclude_lower_left=False, dimness_tolerance=1.0, distance_from_target=900., non_linear_limit=3500, edge_tolerance=125.)

centroided_sources = pat.photometry.centroider(target, sources, restore=False, output_plots=True, gif=False, box_w=5)

pat.photometry.aper_phot(target, centroided_sources, [5], an_in=9)
#pat.photometry.epsf_phot(target, centroided_sources, plots=True)
#pat.photometry.basic_psf_phot(target, centroided_sources, plots=True)

#pat.analysis.lightcurve(target, sources, centroided_sources)
#pat.analysis.speculoos_style_lightcurve(target, sources, output_plots=True, phot_type='aper')
pat.analysis.pines_alc(target, phot_type='aper')