import pines_analysis_toolkit as pat
import pdb

target = '2MASS J17281134+0839590'

sources = pat.photometry.ref_star_chooser(target, restore=True)
pat.analysis.lightcurve(target, sources)

pdb.set_trace()