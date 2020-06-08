import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder
import pdb
import os
import pickle
import glob 
from astropy.visualization import simple_norm
from scipy.spatial import distance
from astropy.modeling import models, fitting
from scipy import signal 
import warnings
import natsort

""" Author: Patrick Tamburo, BU, Dec. 2019
    Purpose: Determines the shift between images and a master image. identify_target_and_reference must have been run first! 
"""

def deshifter(target_name):
    
    def mimir_source_finder(image,sigma_above_bg,fwhm,bpm):
        #Find sources in Mimir images. 

        np.seterr(all='ignore') #Ignore invalids (i.e. divide by zeros)

        #Find stars in the master image.
        daofind = DAOStarFinder(fwhm=fwhm, threshold=sigma_above_bg*stddev, ratio=0.9, exclude_border=True, sky=med)
        new_sources = daofind(image,mask=bpm)

        if new_sources is None:
            print('No sources found! Clouds?')
            pdb.set_trace()

        x_centroids = new_sources['xcentroid']
        y_centroids = new_sources['ycentroid']
        sharpness = new_sources['sharpness']
        fluxes = new_sources['flux']

        #Cut sources that are found within 20 pix of the edges.
        use_pos = np.where((x_centroids > 20) & (x_centroids < 1004) & (y_centroids  > 20) & (y_centroids  < 1004))[0]
        x_centroids = x_centroids[use_pos]
        y_centroids = y_centroids[use_pos]
        sharpness = sharpness[use_pos]
        fluxes = fluxes[use_pos]

        #Also cut on sharpness, this seems to eliminate a lot of false detections.
        use_sharp = np.where(sharpness > 0.4)[0]
        x_centroids  = x_centroids [use_sharp]
        y_centroids  = y_centroids [use_sharp]
        sharpness = sharpness[use_sharp]
        fluxes = fluxes[use_sharp]

        #Finally, cut targets whose y centroids are near y = 512. These are usually bad.
        use_512 = np.where(np.logical_or((y_centroids < 510),(y_centroids > 514)))[0]
        x_centroids  = x_centroids [use_512]
        y_centroids  = y_centroids [use_512]
        sharpness = sharpness[use_512]
        fluxes = fluxes[use_512]
        
        return(x_centroids,y_centroids,fluxes)

    def synthetic_image_maker(x_centroids,y_centroids):
        
        #Construct synthetic images from centroid/flux data.
        synthetic_image = np.zeros((1024,1024))
        sigma = synthetic_fwhm/fwhm_to_sigma

        for i in range(len(x_centroids)):
            #Cut out little boxes around each source and add in Gaussian representations. This saves time. 
            int_centroid_x = int(np.round(x_centroids[i]))
            int_centroid_y = int(np.round(y_centroids[i]))
            y_cut, x_cut = np.mgrid[int_centroid_y-10:int_centroid_y+10,int_centroid_x-10:int_centroid_x+10]
            dist = np.sqrt((x_cut-x_centroids[i])**2+(y_cut-y_centroids[i])**2)
            synthetic_image[y_cut,x_cut] += np.exp(-((dist)**2/(2*sigma**2)+((dist)**2/(2*sigma**2))))
        return(synthetic_image)

    def cross_correlator(check_image,master_synthetic_image):
        #This function finds sources in the current image, makes a synthetic image using those sources,
        #   cross-correlates with the master image, and determines the shift. It executes a while loop,
        #   changing the fwhm to try and detect sources. 
        SEx_fwhm = 7
        while SEx_fwhm < 16:
            (check_x_centroids,check_y_centroids,check_fluxes) = mimir_source_finder(check_image,SEx_sigma,SEx_fwhm,bpm)

            #Create synthetic images from those sources. These functions are defined above.
            check_synthetic_image = synthetic_image_maker(check_x_centroids,check_y_centroids)
            
            #FFT convolve the master and check synthetic images. 
            corr = signal.fftconvolve(master_synthetic_image,check_synthetic_image[::-1,::-1], mode='same')
            
            #Do a 2D Gaussian fit near the pixel with the highest correlation. 
            y_max, x_max = np.unravel_index(np.argmax(corr), corr.shape) 

            #If the maximum pixel is found near an edge, things will crash...set it further away.
            #TODO: come up with a non-idiotic solution. 
            if y_max > 1013: y_max = 1013
            if y_max < 10: y_max = 10
            if x_max > 1013: x_max = 1013
            if x_max < 10: x_max = 10

            y, x = np.mgrid[y_max-10:y_max+10, x_max-10:x_max+10]
            corr_cut =  corr[y,x]
            gaussian_init = models.Gaussian2D(np.max(corr_cut),x_max,y_max,synthetic_fwhm/fwhm_to_sigma,synthetic_fwhm/fwhm_to_sigma,0)
            fit_gauss = fitting.LevMarLSQFitter()
            gaussian = fit_gauss(gaussian_init, x, y,corr_cut)
            x_shift = 512 - gaussian.x_mean.value
            y_shift = 512 - gaussian.y_mean.value

            # If the observing was done properly, we should not measure a shift larger than about 30 pixels
                    #   in x or y. 
            if abs(x_shift) < shift_tolerance and abs (y_shift) < shift_tolerance: 
                break
            else:
                SEx_fwhm = SEx_fwhm + 1
                print('Trying higher source extraction fwhm: fwhm = ',SEx_fwhm)
        
        if abs(x_shift) > shift_tolerance or abs (y_shift) > shift_tolerance:
            print('Found shift larger than 30 pixels in x/y for all values of SEx_fwhm tested. Inspect manually.')
            shift_flags[i] = 1
            pdb.set_trace()
        else:
            return (x_shift,y_shift)

    warnings.filterwarnings('ignore', category=UserWarning, append=True) #This turns of NaN warnings in sigma_clipped_stats, otherwise we'd get a warning every line. 
    fwhm_to_sigma = 2*np.sqrt(2*np.log(2)) #Constant to convert fwhm to sigma. 
    
    synthetic_fwhm = 8 #The fwhm to use for creating synthetic images (don't think the choice matters much).
    shift_tolerance = 30 #The maximum measured shift allowed before things are flagged. If the observing was done properly, you
        #should  never measure a shift > 30 pixels (really more than like 20). 
    SEx_sigma = 3 #Sigma over background for source detection. 

    local_path = '/Users/tamburo/Documents/Data/PINES/'
    reduced_path = local_path+'Objects/'+target_name+'/domeflat_reduced/'
    reduced_files = natsort.natsorted(glob.glob(reduced_path+'*red.fits'))
    
    #Make sure everything is sorted properly.
    # sorting = [int(reduced_files[i].split('.')[1].split('_')[0]) for i in range(len(reduced_files))]
    # sorting_locs = np.argsort(sorting)
    # reduced_files = np.array(reduced_files)[sorting_locs]
    aper_phot_path = local_path+'Objects/'+target_name+'/aper_phot/'
    master_centroids_path = aper_phot_path+'targ_and_refs/'+target_name+'_sources_initial_locs.p'
    bpm = (1-pickle.load(open('/Users/tamburo/Documents/Data/PINES/Calibrations/Bad Pixel Masks/bpm.p','rb'))).astype('bool') 

    output_path = aper_phot_path+'image_shifts/'
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    #Check for previous output. 
    if os.path.exists(output_path+target_name+'_image_shifts.p'):
        old_output = pickle.load(open(output_path+target_name+'_image_shifts.p','rb'))
        old_x_shift = old_output[0]
        old_y_shift = old_output[1]
        old_file_list = old_output[2]
    else:
        old_file_list = []
    
    #Get list of frames.
    all_frame_list = np.array([])
    for i in range(np.size(reduced_files)):
        file_name = reduced_files[i].split('/')[-1]
        all_frame_list = np.append(all_frame_list, file_name)

    #Declare some arrays for saving x/y shifts and flagging bad shifts. 
    save_x_shift = np.zeros(np.size(reduced_files))
    save_y_shift = np.zeros(np.size(reduced_files))
    shift_flags = np.zeros(np.size(reduced_files)) 

    #Create a synthetic image for the master frame.
    master_centroids = pickle.load(open(master_centroids_path,'rb'))
    master_x_centroids = master_centroids[:,0]
    master_y_centroids = master_centroids[:,1]
    master_synthetic_image = synthetic_image_maker(master_x_centroids,master_y_centroids)

    for i in range(np.size(reduced_files)):
        filename = reduced_files[i].split('/')[-1]
        if filename not in old_file_list:
            print('Measuring shifts for ',filename,', ',i+1,' of ',np.size(reduced_files))
            frame = fits.open(reduced_files[i])[0].data
            avg, med, stddev = sigma_clipped_stats(frame, sigma=3.0,maxiters=3)    
            (x_shift,y_shift) = cross_correlator(frame,master_synthetic_image)
            print('X, Y Shift  (pixels): ',np.round(x_shift,1),np.round(y_shift,1))
            save_x_shift[i] = x_shift
            save_y_shift[i] = y_shift
            print('')
        else:
            print(filename,' already has deshifter output, skipping.')
            loc = np.where(np.array(old_file_list) == filename)[0][0]
            save_x_shift[i] = old_x_shift[loc]
            save_y_shift[i] = old_y_shift[loc]
        
    #Save the mesured lists and a list of file names. 
    #
    pickle.dump((save_x_shift,save_y_shift,all_frame_list),open(output_path+target_name+'_image_shifts.p','wb'))
    print('Done deshifting!')
    print('')
