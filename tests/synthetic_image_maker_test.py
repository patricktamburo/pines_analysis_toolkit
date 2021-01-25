import pines_analysis_toolkit as pat 
from astropy.io import fits 
from pathlib import Path
import matplotlib.pyplot as plt 
import pdb 

image_path = Path('/Users/tamburo/Documents/PINES_analysis_toolkit/Objects/2MASS 0026-0943/reduced/20201002.567_red.fits')
image = fits.open(image_path)[0].data
sources = pat.photometry.detect_sources(image_path, seeing_fwhm=2.5, edge_tolerance=90)
source_x = sources['xcenter']
source_y = sources['ycenter']

pat.utils.quick_plot(image)
plt.plot(source_x, source_y, 'rx')
pdb.set_trace()