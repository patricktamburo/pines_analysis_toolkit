import pines_analysis_toolkit as pat
import pdb

targets = ['2MASS J00144919-0838207']
exptimes = [60.]
band = 'J'
flat_date = '20201003'
dark_date = '20201003'
sftp = pat.utils.pines_login()

#Do calibration stuff. 
#for exptime in exptimes:
#    pat.data.dark(dark_date, exptime, sftp=sftp, upload=True, delete_raw=True)
# #pat.data.dome_flat_field(flat_date, band, sftp=sftp, upload=True, delete_raw=True)

#for exptime in exptimes:
#   pat.data.variable_pixels(dark_date, exptime, sftp=sftp, upload=True)
#   pat.data.hot_pixels(dark_date, exptime, sftp=sftp, upload=True)

# #pat.data.dead_pixels(flat_date, band, sftp=sftp, upload=True)

#for exptime in exptimes:
#   pat.data.bpm_maker(flat_date, exptime, band, sftp=sftp, upload=True)

for target in targets:
    pat.data.get_raw_science_files(sftp, target)
    pat.data.reduce(target, sftp=sftp, delete_raw=True, delete_reduced=False, upload=True)
sftp.close()


pdb.set_trace()

sources = pat.photometry.ref_star_chooser(target, restore=True, exclude_lower_left=False, dimness_tolerance=1.2, distance_from_target=300.)
centroided_sources = pat.photometry.centroider(target, sources, restore=True, output_plots=False, box_w=5)

pat.photometry.aper_phot(target, centroided_sources, [4,5,6,7])
#pat.photometry.epsf_phot(target, centroided_sources, plots=True)
#pat.photometry.basic_psf_phot(target, centroided_sources, plots=True)
pat.analysis.lightcurve(target, sources, centroided_sources, phot_type='aper')