import pines_analysis_toolkit as pat 
from astropy.io import fits 
from astropy.wcs import WCS
import pdb 
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from astropy.visualization import ImageNormalize, ZScaleInterval
import matplotlib.pyplot as plt 
plt.ion() 

path = '/Users/tamburo/Documents/PINES_analysis_toolkit/Objects/2MASS 0913+1841/reduced/20200128.002_red_wcs.fits'
hdu = fits.open(path)[0] 
wcs = WCS(hdu.header)

image = interpolate_replace_nans(hdu.data, kernel=Gaussian2DKernel(x_stddev=0.25))
norm = ImageNormalize(image, interval=ZScaleInterval())

fig, ax = plt.subplots(1,1,figsize=(10,10))
ax.imshow(image, origin='lower', norm=norm)
ax.set_xlim(780,920)
ax.set_ylim(730,870)

central_x = 846.74
central_y = 797.90

ax.plot(central_x, central_y, 'kx', ms=10, mew=2)
ax.text(782, 780, '(x , y)      = ({:1.2f},{:1.2f})'.format(central_x,central_y), fontsize=20, color='k', weight='bold')
wcs_pos = wcs.pixel_to_world(central_x, central_y)
ra = wcs_pos.ra.value
dec = wcs_pos.dec.value

ax.text(782, 770, '(RA, Dec) = ({:1.4f},{:1.4f})'.format(ra, dec), fontsize=20, color='k', weight='bold')


pdb.set_trace()