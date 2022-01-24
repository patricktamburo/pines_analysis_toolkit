Useful Functions
================
There are a few functions in the utils module that are regularly helpful in analyzing PINES data, and I'll describe them here. 

jd_utc_to_bjd_tdb()
*******************
The times in Mimir fits headers are reported as Julian dates (JD) in the Coordinated Universal Time (UTC) scale. However, in exoplanet astronomy it has become standard to report times as barycentric Julian dates (BJD) in the 'Temps Dynamique Barycentrique' (TDB) scale, for reasons that are briefly discussed `here <https://astroutils.astronomy.osu.edu/time/bjd_explanation.html>`_ (and in more detail, in Jason Eastman's 2010 paper `here <https://ui.adsabs.harvard.edu/abs/2010PASP..122..935E/abstract>`_ ). This function calculates a BJD_TDB given a JD_UTC. 

pines_log_reader()
******************
This function will read in any PINES csv file into a pandas dataframe, whether an observing log, source detection file, photometry file, etc. 

quick_plot()
************
I use quick_plot constantly. It generates a plot of a 2D array of data with a ZScaleInterval applied (like you can do in DS9). We could call it like this, assuming we're in a directory with a fits file called '20200604.042_red.fits': 

.. code-block:: python 

   from pines_analysis_toolkit.utils import quick_plot as qp 
   from astropy.io import fits

   image = fits.open('20200604.042_red.fits')[0].data
   qp(image)

This generates the following plot: 

.. image:: ../images/qp1.png

In this image, NaNs are indicated as white pixels. We can interpolate these pixels to generate a clearer view of the field: 

.. code-block:: python 

   qp(image, interp=True)

.. image:: ../images/qp2.png

Note that quick_plot() can be run on any 2D array of image data. For example, let's make a cutout around the target and quick plot it: 

.. code-block:: python 
   target_x = 698.6
   target_y = 387.5
   box_w = 30
   cutout = image[int(target_y-box_w/2):int(target_y+box_w/2), int(target_x-box_w/2):int(target_x+box_w/2)]
   qp(cutout, interp=True)

.. image:: ../images/qp3.png