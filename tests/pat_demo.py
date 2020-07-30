import pines_analysis_toolkit as pat

#Set a target.
target = '2MASS J16291840+0335371'

#Set up an sftp connection to the PINES server. 
sftp = pat.utils.pines_login()

#Get reduced data for the target from PINES server.
pat.data.get_reduced_science_files(sftp, target)

#Close the sftp connection. 
sftp.close()

#Identify the target and suitable reference stars using the target's master image. 
sources = pat.photometry.ref_star_chooser(target, restore=False, non_linear_limit=3500, dimness_tolerance=0.9, distance_from_target=800, exclude_lower_left=True)

#Centroid the target and suitable reference stars in every image. 
centroided_sources = pat.photometry.centroider(target, sources, restore=False, plots=False)

#Perform photometry using the centroided sources with 3-, 4-, and 5-pixel radius apertures.
pat.photometry.aper_phot(target, centroided_sources, [3., 4., 5.])

#Make lightcurves. 
pat.analysis.lightcurve(target, sources)

#Look at a movie of the target's centroid position to see if it was properly identified in each image. 
pat.analysis.movie(target, 5)