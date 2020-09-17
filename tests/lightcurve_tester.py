import pines_analysis_toolkit as pat
target = '2MASS J23255604-0259508'
sources = pat.photometry.ref_star_chooser(target, restore=True)
pat.analysis.lightcurve(target, sources)