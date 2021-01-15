import pines_analysis_toolkit as pat
import pdb

targets = ['2MASS J00191165+0030176']
#targets = ['SIMP 0136'] #Declare a target. This is SIMP J013656.5+093347.
exptimes = [10.] #Exposure time of the images in seconds. 
flat_date = '20201003' #The date that dome flat calibration data were taken for the target. 
dark_date = '20201003' #The date that dark calibration data were taken for the target.
band = 'J'

#sftp = pat.utils.pines_login()

#for exptime in exptimes:
#     pat.data.dark(dark_date, exptime, upload=True, delete_raw=True, sftp=sftp)
#pat.data.dome_flat_field(flat_date, band, sftp=sftp, upload=True, delete_raw=True)

#for exptime in exptimes:
#    pat.data.variable_pixels(dark_date, exptime, upload=True, sftp=sftp)
#    pat.data.hot_pixels(dark_date, exptime, box_l=5, upload=True, sftp=sftp)

#pat.data.dead_pixels(flat_date, band, upload=False, sftp=sftp)

#for exptime in exptimes:
#    pat.data.bpm_maker(flat_date, dark_date, exptime, band, upload=True, sftp=sftp)

#for target in targets:
#    pat.data.get_raw_science_files(sftp, target)
#    pat.data.reduce(target, delete_raw=True, delete_reduced=False, upload=True, sftp=sftp)

#sftp.close()


target = targets[0]
sources = pat.photometry.ref_star_chooser(target, restore=True, source_detect_image_ind=120, 
            exclude_lower_left=False, dimness_tolerance=0.1, distance_from_target=900., non_linear_limit=3300, 
            edge_tolerance=30., source_detect_plot=False)


centroided_sources = pat.photometry.centroider(target, sources, restore=True, output_plots=True, gif=False, box_w=7)

#pat.photometry.aper_phot(target, centroided_sources, [6], an_in=9)
#pat.photometry.epsf_phot(target, centroided_sources, plots=True)
#pat.photometry.basic_psf_phot(target, centroided_sources, plots=True)

#pat.analysis.lightcurve(target, sources, centroided_sources)
#pat.analysis.speculoos_style_lightcurve(target, sources, output_plots=True, phot_type='aper')
pat.analysis.pines_alc(target, phot_type='aper', mode='night')