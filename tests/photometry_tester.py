import pines_analysis_toolkit as pat

#Declare a target.
target = '2MASS J16291840+0335371'

#Set up connection to pines.bu.edu.
sftp = pat.utils.pines_login()

#Grab reduced data for target. 
pat.data.get_reduced_science_files(sftp, target)

#Close sftp connection. 
sftp.close()

#Detect target and suitable reference stars. 
targ_and_refs = pat.photometry.ref_star_chooser(target, restore=False, lower_left=True)

#Centroid target and references. 
centroided_sources = pat.photometry.centroider(target, targ_and_refs, restore=False)

#Do aperture photometry on target and references. 
pat.photometry.aper_phot(target, centroided_sources, [5.,6.,7.,8.,9.])
