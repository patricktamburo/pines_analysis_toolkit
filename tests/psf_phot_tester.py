import pines_analysis_toolkit as pat

target = '2MASS J23255604-0259508'

#Detect target and suitable reference stars. 
targ_and_refs = pat.photometry.ref_star_chooser(target, source_detect_image='', restore=True, exclude_lower_left=True, dimness_tolerance=1.0, distance_from_target=800., non_linear_limit=3500)

#Centroid target and references. 
centroided_sources = pat.photometry.centroider(target, targ_and_refs, restore=True, plots=True)

#Do aperture photometry on target and references. 
pat.photometry.psf_phot(target, centroided_sources, plots=False)