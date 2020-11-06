import pines_analysis_toolkit as pat
import pdb

target = '2MASS J01365662+0933473'
exptime = 30.
band = 'J'
flat_date = '20201102'
dark_date = '20201102'

#sftp = pat.utils.pines_login()

#Do calibration stuff. 

#for exp in exptime:
#    pat.data.dark(dark_date, exp, sftp=sftp, upload=True, delete_raw=False)

#pat.data.dome_flat_field(flat_date, band, sftp=sftp, upload=True, delete_raw=False)


# for exp in exptime:
#     pat.data.variable_pixels(dark_date, exp, sftp=sftp, upload=True)
#     pat.data.hot_pixels(dark_date, exp, sftp=sftp, upload=True)

# pat.data.dead_pixels(flat_date, band, sftp=sftp, upload=True)

# for exp in exptime:
#     pat.data.bpm_maker(flat_date, exp, band, sftp=sftp, upload=True)

# sftp.close()

# pdb.set_trace()

#pat.data.get_raw_science_files(sftp, target)

pat.data.reduce(target, upload=False, delete_raw=False, delete_reduced=False)

#pdb.set_trace()
sources = pat.photometry.ref_star_chooser(target, restore=True, source_detect_image='20201104.009_red.fits', exclude_lower_left=False, dimness_tolerance=0.4, distance_from_target=1000.)
centroided_sources = pat.photometry.centroider(target, sources, restore=False, output_plots=False, box_w=17)

pat.photometry.aper_phot(target, centroided_sources, [5,6,7,8], an_in=10, an_out=25)
#pat.photometry.epsf_phot(target, centroided_sources, plots=True)
#pat.photometry.basic_psf_phot(target, centroided_sources, plots=True)
#pat.analysis.speculoos_style_lightcurve(target, phot_type='aper')
#pdb.set_trace()
pat.analysis.lightcurve(target, sources, centroided_sources, phot_type='aper')