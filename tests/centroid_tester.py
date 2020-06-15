import pines_analysis_toolkit as pat
import pandas as pd 
import pdb

target = '2MASS J16291840+0335371'

targ_and_refs = pd.read_csv('/Users/tamburo/Documents/PINES_analysis_toolkit/Objects/2MASS 1629+0335/sources/target_and_references.csv')
centroided_sources = pat.photometry.centroider(target, targ_and_refs)
pdb.set_trace()