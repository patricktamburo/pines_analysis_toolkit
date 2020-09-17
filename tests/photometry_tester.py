import pines_analysis_toolkit as pat
import matplotlib.pyplot as plt
import pdb

target = '2MASS J23255604-0259508'

# #Set up connection to pines.bu.edu.
# sftp = pat.utils.pines_login()

# #Grab raw data for target. 
# pat.data.get_raw_science_files(sftp, target)

# #Reduce data
# pat.data.reduce(target, delete_raw=False, upload=False) 

# #Close sftp connection. 
# sftp.close()

#Detect target and suitable reference stars. 
targ_and_refs = pat.photometry.ref_star_chooser(target, source_detect_image='', restore=False, exclude_lower_left=True, dimness_tolerance=1.0, distance_from_target=800., non_linear_limit=3500)

#Centroid target and references. 
centroided_sources = pat.photometry.centroider(target, targ_and_refs, restore=False, plots=True)

#Do aperture photometry on target and references. 
pat.photometry.aper_phot(target, centroided_sources, [3,4,5,6,7], plots=False)
