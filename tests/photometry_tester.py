import pines_analysis_toolkit as pat
import pdb

target = '2MASS J16291840+0335371'

targ_and_refs = pat.photometry.ref_star_chooser(target, restore=True)
centroided_sources = pat.photometry.centroider(target, targ_and_refs, restore=True)
aper_phot = pat.photometry.aper_phot(target, centroided_sources, [5.,6.,7.,8.,9.])

pdb.set_trace()