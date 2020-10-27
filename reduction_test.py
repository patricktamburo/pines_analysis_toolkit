import pines_analysis_toolkit as pat
import pdb

target = '2MASS J04384965+4349164'
exptime = 15.
band = 'J'
flat_date = '20201003'
dark_date = '20200129'
#sftp = pat.utils.pines_login()

#Do calibration stuff. 
#pat.data.dark(dark_date, exptime, dark_start=1, dark_stop=20)
#pat.data.dome_flat_field(flat_date, band)
#pat.data.variable_pixels(dark_date, exptime)
#pat.data.hot_pixels(dark_date, exptime)
#pat.data.dead_pixels(flat_date, band)
#pat.data.bpm_maker(flat_date, exptime, band)
#pat.data.get_raw_science_files(sftp, target)
#sftp.close()

#pat.data.reduce(target, upload=False, delete_raw=False, delete_reduced=False)

#pdb.set_trace()

sources = pat.photometry.ref_star_chooser(target, restore=True, exclude_lower_left=False, dimness_tolerance=1.2, distance_from_target=300.)
centroided_sources = pat.photometry.centroider(target, sources, restore=True, plots=True, box_w=5)

pat.photometry.aper_phot(target, centroided_sources, [4,5,6,7])
#pat.photometry.epsf_phot(target, centroided_sources, plots=True)
#pat.photometry.basic_psf_phot(target, centroided_sources, plots=True)
pat.analysis.lightcurve(target, sources, phot_type='aper')