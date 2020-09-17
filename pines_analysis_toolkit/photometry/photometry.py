import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder
from photutils import find_peaks
import pdb
import os
import pickle
import glob 
from astropy.visualization import simple_norm
from scipy.spatial import distance
import matplotlib
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import aperture_photometry
from scipy.optimize import curve_fit
from astropy.time import Time
import natsort

#----------------------------------------------------------------USER INPUTS---------------------------------------------------------------------------------
def photometry(target_name,all_apertures,an_in=12.,an_out=30.):
    """ 
        AUTHOR: 
            Patrick Tamburo, BU, Dec. 2019
        PURPOSE: 
            Does aperture photometry on Mimir images
        INPUTS:
            target_name (str): The name of your target, which must match the data analysis folder name. 
            local_path (str): The analysis directory for this target.
            all_apertures (list): Radii of photometric aperetures to use. An output file will be made for each.
            an_in (float): The inner radius annulus, in pixels. By default, 12.
            an_out (float): The outer radius annulus, in pixels. By default, 30.
        OUTPUTS:
    """
    
    def phot(position):
        #Define the aperture. 
        phot_aperture = CircularAperture(position,r=ap)
        #Define the annulus.
        phot_annulus = CircularAnnulus(position, r_in=an_in, r_out=an_out)

        raw_ap_flux = aperture_photometry(frame,phot_aperture)
        raw_an_flux = aperture_photometry(frame,phot_annulus)

        mask = phot_annulus.to_mask()
        # Construct a histogram of pixel values contained entirely within
        # the annulus.
        for ii in range(0, len(mask.data)):
            for jj in range(0, len(mask.data[0])):
                if mask.data[ii, jj] < 1.0:
                    mask.data[ii, jj] = 0.
        data_cutout = mask.cutout(frame) * mask


        bg_pix_vals = []
        # Take the non-zero values into a list.data_cutout = mask[0].cutout(data) * mask_obj
        for ii in range(0, len(data_cutout)):
            for jj in range(0, len(data_cutout[0])):
                if data_cutout[ii, jj] != 0 and np.isnan(data_cutout[ii,jj]) == False:
                    bg_pix_vals.append(data_cutout[ii, jj])

        bg_pix_vals = np.array(bg_pix_vals)
        median_bg = np.median(bg_pix_vals)
        std_bg = np.std(bg_pix_vals)
        # Remove >3-sigma outliers (cosmic ray hits, bad pixels, contamination, etc.)
        outlier_locs = np.where(bg_pix_vals > median_bg + 3. * std_bg)[0]
        bg_filtered_vals = np.delete(bg_pix_vals, outlier_locs)

        #Remove negative pixels
        outlier_locs = np.where(bg_filtered_vals < 0)[0]
        bg_filtered_vals = np.delete(bg_filtered_vals, outlier_locs)

        an_mean = np.median(bg_filtered_vals)

        #Now, correct the raw aperture flux
        ap_flux = raw_ap_flux['aperture_sum'] - an_mean*phot_aperture.area

        return(ap_flux, an_mean)

    all_apertures = np.array(all_apertures)

    local_path = '/Users/tamburo/Documents/Data/PINES/'
    reduced_path = local_path+'Objects/'+target_name+'/domeflat_reduced/'
    reduced_files = natsort.natsorted(glob.glob(reduced_path+'*red.fits'))

    #Make sure everything is sorted properly.
    # sorting = [int(reduced_files[i].split('.')[0].split('/')[-1]+reduced_files[i].split('.')[1].split('_')[0]) for i in range(len(reduced_files))]
    # sorting_locs = np.argsort(sorting)
    # reduced_files = np.array(reduced_files)[sorting_locs]

    aper_phot_path = local_path+'Objects/'+target_name+'/aper_phot/'
    target_and_ref_path = aper_phot_path+'targ_and_refs/'
    image_shift_path = aper_phot_path+'image_shifts/'
    centroid_path = aper_phot_path+'centroids/'
    bpm = (1-pickle.load(open('/Users/tamburo/Documents/Data/PINES/Calibrations/Bad Pixel Masks/bpm.p','rb'))).astype('bool') 

    output_path = aper_phot_path+'photometry/'
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    #Get list of files
    all_frame_list = np.array([])
    for i in range(np.size(reduced_files)):
        file_name = reduced_files[i].split('/')[-1]
        all_frame_list = np.append(all_frame_list, file_name)
    num_files = len(reduced_files)

    #Get image positions of target and references.
    positions = pickle.load(open(centroid_path+target_name+'_centroids.p','rb'))
    x_positions = positions['X pos']
    y_positions =  positions['Y pos']
    x_subframe_positions = positions['X subframe pos']
    y_subframe_positions = positions['Y subframe pos']
    num_targets = len(x_positions) 


    #Get seeing values/centroid positions for saving purposes.
    seeing = pickle.load(open(centroid_path+target_name+'_seeing.p','rb'))
    x_seeing = seeing['X seeing']
    y_seeing = seeing['Y seeing']

    for ap in all_apertures:

        #Check for existing output.
        output_filename = output_path+target_name+'_ap='+str(ap)+'_photometry.p'
        if os.path.exists(output_filename):
            old_output = pickle.load(open(output_path+target_name+'_ap='+str(ap)+'_photometry.p','rb'))
            old_time = old_output['Time']
            old_flux = old_output['Flux']
            old_airmass =old_output['Airmass']
            old_background = old_output['Background']
            old_file_list = old_output['Old frame list']
        else:
            old_file_list = []
        
        #Create some empty arrays for outputting data. 
        save_airmass = np.zeros(num_files)
        save_flux = np.zeros((num_targets,num_files)) #The measured flux for the target/references
        save_bkg = np.zeros((num_targets,num_files))
        save_time = np.zeros(num_files) #The time of each frame 
        
        #Loop over all files. 
        for i in range(num_files):
            filename = reduced_files[i].split('/')[-1]
            if filename not in old_file_list:
                print('Doing ap = ',ap,' photometry for',filename,', ',i,' of ',np.size(reduced_files))
                frame = fits.open(reduced_files[i])[0].data
                header = fits.open(reduced_files[i])[0].header
                save_time[i] = Time(header['DATE-OBS']).jd 
                save_airmass[i] = header['Airmass']

                #Loop over all targets in image.
                for j in range(num_targets):
                    position = (x_positions[j,i],y_positions[j,i])
                    (save_flux[j,i], save_bkg[j,i]) = phot(position)

                    pdb.set_trace()

                #End loop over targets.
            else:
                #If the photometry was done previously, just read in the old values. 
                print('Ap = ',str(ap),' photometry already exists for',filename,', skipping.')
                loc = np.where(np.array(old_file_list) == filename)[0][0]
                save_time[i] = old_time[loc]
                save_airmass[i] = old_airmass[loc]
                save_bkg[:,i] = old_background[:,loc]
                save_flux[:,i] = old_flux[:,loc]
            #End loop over files.


        output_dict = {'Time':save_time,'Flux':save_flux,'Airmass':save_airmass,'Background':save_bkg,
                        'X position':x_positions,'Y position':y_positions,
                        'X subframe position':x_subframe_positions,'Y subframe position':y_subframe_positions,
                        'X seeing':x_seeing,'Y seeing':y_seeing,'Old frame list':all_frame_list}
        pickle.dump(output_dict,open(output_filename,'wb'))
        print('Done photometry for ap = ',ap,'!')
        print('')
    #End loop over apertures.
