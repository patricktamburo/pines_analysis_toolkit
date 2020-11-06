import pines_analysis_toolkit as pat
import pysftp
import pdb 

target = '2MASS J23062928-0502285'
exptime = 10.
band = 'J'
calibration_date = '20200907'

sftp = pat.utils.pines_login()

#Do calibration stuff. 
pat.data.dark(calibration_date, exptime, upload=True, sftp=sftp)
#pat.data.dome_flat_field(calibration_date, band, upload=True, sftp=sftp)
pat.data.variable_pixels(calibration_date, exptime, upload=True, sftp=sftp)
pat.data.hot_pixels(calibration_date, exptime, upload=True, sftp=sftp)
#pat.data.dead_pixels(calibration_date, band, upload=True, sftp=sftp)
pat.data.bpm_maker(calibration_date, exptime, band, upload=True, sftp=sftp)

sftp.close()

pdb.set_trace()

#Reduce the data.
#pat.data.reduce(target)

# sources = pat.photometry.ref_star_chooser(target, restore=True, source_detect_image='20200905.197_red.fits', dimness_tolerance=0.25, exclude_lower_left=True, distance_from_target=800)

# centroided_sources = pat.photometry.centroider(target, sources, restore=True, plots=True)

# pat.photometry.epsf_phot(target, centroided_sources, plots=True)

# pat.analysis.lightcurve(target, sources, phot_type='psf')