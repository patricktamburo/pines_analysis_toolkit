import pines_analysis_toolkit as pat 

target = '2MASS J01365662+0933473'
sources = pat.photometry.ref_star_chooser(target, dimness_tolerance=0.1, restore=True)
centroided_sources = pat.photometry.centroider(target, sources, restore=True)
pat.photometry.aper_phot(target, centroided_sources, [4,5])
pat.analysis.speculoos_style_lightcurve(target)