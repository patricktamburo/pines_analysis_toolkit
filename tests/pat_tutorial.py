import pines_analysis_toolkit as pat
import pdb

'''
-------------------------------INSTRUCTIONS BEFORE YOU START------------------------------------------------
IF AN EARLIER VERSION OF PINES_ANALYSIS_TOOLKIT IS INSTALLED, YOU MUST REMOVE IT AND REINSTALL TO SEE UPDATES!

Do this by opening a terminal and typing the following:
   pip3 uninstall pines_analysis_toolkit
...and hit y when prompted. (Your command for this may use pip instead of pip3, but pip3 is what works on my computer)

Now, pull the new version of the repository to wherever you have the pines_analysis_toolkit files downloaded, by navigavting there with a terminal and typing:
    git pull origin master

This assumes you have already set up the Github repository as your 'origin'; if you haven't, first do that by typing: 
    git remote add origin https://github.com/patricktamburo/pines_analysis_toolkit.git

Now, REINSTALL the pines_analysis_toolkit by typing:
    python3 setup.py install

If it installs successfully, then you're able to run this code to test out some of the features. 
-------------------------------------------------------------------------------------------------------------
'''

#Declare a target that you want reduced science data for, using its full 2MASS name.
target = '2MASS J16291840+0335371'
pat.data.get_reduced_science_files(target)

#Find the initial location for the target in the reduced images, and get a set of suitable reference stars. 
#The Mimir "bars" issue is present in some of the 2M 1629+0335 data, so I turn on the lower_left flag to automatically discard references in the lower left quadrant.
sources = pat.photometry.ref_star_chooser(target, lower_left=True)

#Get positions of these sources in all of the images. 
centroided_sources = pat.photometry.centroider(target, sources)


pdb.set_trace()