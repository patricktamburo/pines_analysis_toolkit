import PINES as p
import pdb
import os
import time

def photometry_main(target_name,aperture_radii,guess_position=(705,386),source_frame_id=0,run_identify_sources=0,run_deshifter=1,run_centroider=1,run_photometry=1):
    #Get the target and reference star initial positions.
    if run_identify_sources:
        print('Running identify_sources.')
        print('')
        time.sleep(1)
        p.photometry.identify_sources(target_name,guess_position,source_frame_id=50,SEx_fwhm=8)
    else:
        print('Restoring identify_sources output.')
        print('')
        time.sleep(1)

    #Measure shifts between the master image and other images. 
    if run_deshifter:
        print('Running deshifter.')
        print('')
        time.sleep(1)
        p.photometry.deshifter(target_name)
    else:
        print('Restoring deshifter output.')
        print('')
        time.sleep(1)

    #Get centroid positions of target and references. 
    if run_centroider:
        print('Running centroider.')
        print('')
        time.sleep(1)
        p.photometry.centroider(target_name)
    else:
        print('Restoring centroider output.')
        print('')
        time.sleep(1)

    #Get centroid positions of target and references. 
    if run_photometry:
        print('Running photometry.')
        print('')
        time.sleep(1)
        p.photometry.photometry(target_name,aperture_radii)

