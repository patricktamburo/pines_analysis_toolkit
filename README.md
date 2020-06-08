# pines_analysis_toolkit

[PINES website](https://pines.bu.edu/)

This is a package for working with PINES data.

INSTALLATION
  1. Install git on your computer. 
  2. Create a directory where you want pines_analysis_toolkit to be installed. When writing new features for the toolkit, you will    do so in this directory.
  3. Navigate a terminal to that directory, and type the following: 
  
        git init
        
        git remote add origin https://github.com/patricktamburo/pines_analysis_toolkit.git
        
        git pull origin master
        
  4. Check to make sure files have appeared in the directory.
  5. Install the package locally by typing the following:
      python setup.py install
    OR
      python3 setup.py install
      ...whichever is the appropriate command for python 3 on your computer. 
  6. Check your Documents folder to see if the PINES_analysis_toolkit directory was created. All PINES data that you analyze will be placed here. 
  7. Start python, and type import pines_analysis_package. 
  8. If it imports without throwing any errors, you're done!

MODULES (all still under development)
 
    data - programs to grab raw data, science data, calibrations, logs, etc. 
  
    reduction - programs for reducing raw data
  
    photometry - programs for performing photometry on reduced data
  
    analysis - programs for making lightcurves, analyzing variability, etc.
    
    utils - programs for various quality-of-life things
