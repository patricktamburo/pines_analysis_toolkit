import pines_analysis_toolkit as pat
import pdb 

target = '2MASS J23255604-0259508'

sources = pat.photometry.ref_star_chooser(target, restore=True)
centroided_sources = pat.photometry.centroider(target, sources, restore=True)
pat.data.image_cleaner(target, centroided_sources)
pdb.set_trace()   