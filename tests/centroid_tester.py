import pines_analysis_toolkit as pat
import pandas as pd 
import pdb

target = '2MASS J16291840+0335371'

targ_and_refs = pat.photometry.ref_star_chooser(target, restore=True)
centroided_sources = pat.photometry.centroider(target, targ_and_refs, restore=True)
pdb.set_trace()