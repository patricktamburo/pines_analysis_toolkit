Introduction and Installation
============
PINES analysis toolkit, or PAT, is a python package for interacting with data from the Perkins INfrared Satellite Survey, a search for transiting exosatellites around a sample of almost 400 L and T dwarfs. It contains all of the tools necessary to go from raw images to final light curves. The pipeline is still under active development, so check the project's `Github <https://github.com/patricktamburo/pines_analysis_toolkit>`_ frequently!

Installation
************

PINES_analysis_toolkit Directory Structure
******************************************
When PAT is installed, it will by default create a folder in your ~/Documents/ folder called 'PINES_analysis_toolkit'. Generally, we refer to this folder as your 'PINES path'. Everything related to PINES data will reside in your PINES path. It has the following structure:

::

    PINES_analysis_toolkit
    ├── Calibrations
    │  ├── Bad Pixel Masks
    │  ├── Darks
    │  └── ...
    ├── Logs
    ├── Misc          
    └── Objects          
        ├── 2MASS 0014-0838
          ├── 2mass0014-0838.profile
          ├── analysis
          ├── aper_phot
          ├── output
          ├── psf_phot
          ├── raw
          ├── reduced
          └── sources
        ├── 2MASS 0019+0030
        └── ...

As the names suggest, calibration files will be created in Calibrations, observing logs will be saved to Logs, and miscellaneous supporting documents are in Misc. Most of the time, however, we'll be worrying about 'Objects' directory, which contains subdirectories for every object that you have done PAT analysis on. When you first install PAT, nothing will exist in your Objects folder. We'll work on populating it in the following Tutorial.