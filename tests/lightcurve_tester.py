import pines_analysis_toolkit as pat
import pdb

target = '2MASS J16291840+0335371'
sources = pat.photometry.ref_star_chooser(target, restore=True)
lc = pat.analysis.lightcurve(target, sources)

pdb.set_trace()