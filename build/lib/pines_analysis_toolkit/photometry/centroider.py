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
from scipy.optimize import curve_fit
import natsort

""" Author: Patrick Tamburo, BU, Dec. 2019
    Purpose: Determines centroid locations for your target and references in Mimir images. 
        Also measures x/y seeing in each image.
        identify_target_and_reference.py and deshifter.py must have been run first!

"""


def centroider(target_name,sub_frame_l=40,SEx_fwhm=6.0,SEx_sigma=3.0,plot_subframes=False):
    """ 
        AUTHOR: 
            Patrick Tamburo, BU, Dec. 2019
        PURPOSE: 
            Determines the centroid positions of your target/references (as determined in identify_sources()) using de-shifted 
                subframes. Also measures the seeing in x and y around your target in each image. 
        INPUTS:
            target_name (str): The name of your target, which must match the data analysis folder name. 
            local_path (str): The analysis directory for this target.
            sub_frame_l (int): The size of the subframe (in pixels) used to centroid a target. By default, 40.  
            SEx_fwhm (float): The fwhm used in DAOStarFinder to detect stars. By default, set to 6. 
            SEx_sigma (float): The sigma above background used in DAOStarFinder to detect stars. By default, set to 3. 
            plot_subframes (bool): If true, will plot the subframes/centroids for your targets in each image. Good for debugging!
        OUTPUTS:
            Writes two files to local_path/aper_phot/centroids:
                1) target_name_centroids.p: A pickle file of the centroid positions determined for your target/references in each image.
                2) target_name_seeing.p: A pickle file containing the measured seeing in each image. 
    """
    
    def get_source(frame,bpm,cen_position,id):
        #Run a source extraction on a subframe centered around the position of the target
        # in the source detection frame. The subframe is a box of side length defined here
        # The frame reads the fits file such that the bottom left corner of the file is 
        # index [0,0].
        avg2, med2, stddev2 = sigma_clipped_stats(frame,sigma=3.0,maxiters=3) 
        if plot_subframes:
            #Plot this subframe in the already launched figure. 
            ax1.clear()
            ax1.imshow(frame, origin='lower',vmin=med2,vmax=med2+12*stddev2)
            if id == 0:
                plt.title(filename+'\nTarget, Image '+str(i+1)+' of '+str(num_files))
            else:
                plt.title(filename+'\nRef '+str(id)+', Image '+str(i+1)+' of '+str(num_files))
            plt.ion()
        
        #Identify bad pixels in the sub_frame, ignore them in daostarfinder
        #new_bpm = np.zeros(np.shape(frame)).astype('bool')
        #new_bpm[np.where(frame < med2-7*stddev2)] = True 
        #TODO: change threshold to a variable, instead of locking it to 5. 
        #daofind = DAOStarFinder(fwhm=SEx_fwhm, threshold=5*stddev2,roundlo=-1000,roundhi=1000,sharplo=-1000,sharphi=1000)

        # if (i == 0) or (reduced_files[i].split('/')[-1].split('.')[0] != reduced_files[i-1].split('/')[-1].split('.')[0]):
        #     #Generally worse seeing at the start of the nights.
        #     daofind = DAOStarFinder(fwhm=SEx_fwhm, threshold=5*stddev2,roundlo=-1000,roundhi=1000,sharplo=-1000,sharphi=1000,ratio=0.9,sky=med2)
        # else:
        #     #If you have info from previous images, use that to find sources. 
        #     daofind = DAOStarFinder(fwhm=np.mean([save_x_seeing[i-1],save_y_seeing[i-1]])/plate_scale, threshold=4*stddev2,roundlo=-1000,roundhi=1000,sharplo=-1000,sharphi=1000,ratio=0.9,sky=med2)

        # new_sources = daofind(frame)

        SEx_fwhm = 3
        max_SEx_fwhm = 16
        while SEx_fwhm < max_SEx_fwhm:
            daofind = DAOStarFinder(fwhm=SEx_fwhm, threshold=4*stddev2,roundlo=-1000,roundhi=1000,sharplo=-1000,sharphi=1000,ratio=0.9,sky=med2)
            new_sources = daofind(frame)

            #If no sources are found, then the source must either be unfocused or too weak 
            # to be detected. Or there may have been a huge shift, which is dealt with in 
            # another portion of the code. If no detection is made with the orignial 
            # DAOFIND arguments try several other calls meant to addressed defocused or
            # low signal/high background fields

            ###DON"T NEED TO DO THIS IF ALREADY LOOPING OVER FWHM
            # if new_sources is None:
            #     #First calculate the max of the sub frame to help diagnose the issue
            #     sub_frame_max = np.amax(frame)
            #     #If the max is still well above the median, then the source is likely out
            #     # of focus. So we will run another source detection with a slightly lower
            #     # threshold and a larger FWHM
                
            #     if sub_frame_max > 5.:
            #         daofind = DAOStarFinder(fwhm=1.25*SEx_fwhm,threshold=(0.75*SEx_sigma)*stddev2)
            #         new_sources = daofind(frame - med2)
            #         if new_sources is not None:
            #             print('Source found in DAOStarFinder second attempt.')
            #     #Otherwise, the signal is probably low, so simply reduce the fwhm
            #     else:
            #         daofind = DAOStarFinder(fwhm=0.75*SEx_fwhm, threshold=(0.75*SEx_sigma)*stddev2)
            #         new_sources = daofind(frame - med2)
            #         if new_sources is not None:
            #             print('Source found in DAOStarFinder second attempt.')
                    
            #If still no sources are found, use find_peaks to locate the flux peak above
            # the given threshold. Force a single source to be returned. If this method 
            # does not return a source, then the flux did not surpass the threshold, and
            # the following block will be executed.
            
            # #TODO: findpeaks can return sub-pixel centroids if centroid_func = centroid_2dg is passed. 
            # if new_sources is None:
            #     #Run find_peaks
            #     #new_sources = find_peaks(frame-med2,threshold=SEx_sigma*stddev2,box_size=sub_frame_l,subpixel=True,npeaks=1,mask=new_bpm)
            #     pdb.set_trace()
            #     new_sources = find_peaks(frame-med2,threshold=SEx_sigma*stddev2,box_size=sub_frame_l,subpixel=True,npeaks=1)

            #     if new_sources is not None:
            #         print('Source found in find_peaks first try!')
            #         new_sources.rename_column('x_peak','xcentroid')
            #         new_sources.rename_column('y_peak','ycentroid')
            #         new_sources.rename_column('peak_value','peak')

            # #If still no sources are found, use find_peaks again but on a much lower 
            # # threshold. Again, force only 1 peak to be found
            # if  new_sources is None:
            #     #Run find_peaks
            #     #new_sources = find_peaks(frame-med2,threshold=0.5*SEx_sigma*stddev2,box_size=sub_frame_l,subpixel=True,npeaks=1,mask=new_bpm)
            #     new_sources = find_peaks(frame-med2,threshold=0.5*SEx_sigma*stddev2,box_size=sub_frame_l,subpixel=True,npeaks=1)

            #     if new_sources is not None:
            #         print('Source found in find_peaks second try!')
            #         new_sources.rename_column('x_peak','xcentroid')
            #         new_sources.rename_column('y_peak','ycentroid')
            #         new_sources.rename_column('peak_value','peak')	

            #If *no* sources wre found, return the center of the frame. 
            if new_sources is None:
                print('No sources were found, returning the center of the frame and flagging. Setting seeing to previous value.')
                save_centroiding_flags[id,i] = 1
                save_x_seeing[i] = save_x_seeing[i-1]
                save_y_seeing[i] = save_y_seeing[i-1]
                if SEx_fwhm == max_SEx_fwhm-1:
                    #If it's the last tested value, return the center
                    return ((cen_position[0],cen_position[1]),((sub_frame_l+1)/2,(sub_frame_l+1)/2), frame)
                else:
                    #Else continue on to the next fwhm value
                    SEx_fwhm = SEx_fwhm + 1
                    #pdb.set_trace()
                    continue

                
            #If only one source is found, this is the target. Determine the centroid of this 
            # new source relative to the full frame and then create its offset
            if np.size(new_sources) == 1:
                #For the full frame

                new_targ_x = new_sources[0]['xcentroid'] + ((cen_position[0])-sub_frame_l//2) - (cen_position[0]-int(cen_position[0]))
                new_targ_y = new_sources[0]['ycentroid'] + ((cen_position[1])-sub_frame_l//2) - (cen_position[1]-int(cen_position[1]))

                #Otherwise, find the offset for this new position
                offsets = [new_targ_x-cen_position[0],new_targ_y-cen_position[1]]
                
                #Modified this to catch large centroid deviations, not just those outside of the box. Don't expect the target to deviate by more than +/- 5 pixels
                #	because of the cross-correlation shifting we did. 
                if (new_sources[0]['xcentroid'] > sub_frame_l/2+10) or (new_sources[0]['xcentroid'] < sub_frame_l/2-10) \
                    or (new_sources[0]['ycentroid'] > sub_frame_l/2+10) or (new_sources[0]['ycentroid'] < sub_frame_l/2-10):
                    print('Large deviation measured for target ',id,'in ',filename,', flagging.')
                    save_centroiding_flags[id,i] = 1
                    if SEx_fwhm == max_SEx_fwhm:
                        pdb.set_trace()
                        return ((cen_position[0],cen_position[1]),((sub_frame_l+1)/2,(sub_frame_l+1)/2),frame)
                    else:
                        SEx_fwhm = SEx_fwhm + 1
                        continue
                
                else:
                    if id == 0:
                        #Estimate the seeing around the target 
                        (fwhm_x,fwhm_y) = seeing_estimator(frame,int(round(new_sources[0]['xcentroid'])),int(round(new_sources[0]['ycentroid'])))
                        save_x_seeing[i] = fwhm_x
                        save_y_seeing[i] = fwhm_y
                    #If a shift was applied correctly, should be within 1.5 pixels of the box center. 
                    if (abs(new_sources[0]['xcentroid']-sub_frame_l/2) > 2.5) or ((abs(new_sources[0]['ycentroid']-sub_frame_l/2) > 2.5) ):
                        #print('Non-centered centroid determined for id = ',str(id),', trying higher fwhm.')
                        if (SEx_fwhm == max_SEx_fwhm-1):
                            print('Non-centered centroid determined for target ',id,', flagging.')
                            save_centroiding_flags[id,i] = 1
                            return ((cen_position[0],cen_position[1]),((sub_frame_l+1)/2,(sub_frame_l+1)/2),frame)
                        else:
                            SEx_fwhm = SEx_fwhm + 1
                            continue
                    else:
                        if plot_subframes:
                            #Plot the measured centroid position. 
                            centroid1, = ax1.plot([new_sources[0]['xcentroid']], [new_sources[0]['ycentroid']], 'rx')
                            plt.pause(0.1)
                            centroid1.remove()
                            ax1.clear()
                        #if SEx_fwhm == max_SEx_fwhm:
                        #    pdb.set_trace()
                        return ((new_targ_x,new_targ_y),(new_sources[0]['xcentroid'],new_sources[0]['ycentroid']), frame)

            if np.size(new_sources) > 1:         
                distances = np.sqrt((new_sources['xcentroid']-sub_frame_l//2)**2+(new_sources['ycentroid']-sub_frame_l//2)**2)
                min_distance_index = np.where(distances == min(distances))[0][0]

                new_targ_x = new_sources[min_distance_index]['xcentroid']+ ((cen_position[0])-sub_frame_l//2) - (cen_position[0]-int(cen_position[0]))
                new_targ_y = new_sources[min_distance_index]['ycentroid']+ ((cen_position[1])-sub_frame_l//2) - (cen_position[1]-int(cen_position[1]))

                
                #Otherwise, find the offset for this new position	
                offsets = [new_targ_x-cen_position[0],new_targ_y-cen_position[1]]

                
                

                #if (new_sources[min_distance_index]['xcentroid'] > sub_frame_l) or (new_sources[min_distance_index]['ycentroid'] > sub_frame_l) or (new_sources[min_distance_index]['xcentroid'] < 0) or (new_sources[min_distance_index]['ycentroid'] < 0):
                #Again, modified to catch large deviations, not just those falling out of frame. 
                if (new_sources[min_distance_index]['xcentroid'] > sub_frame_l/2+20) or (new_sources[min_distance_index]['xcentroid'] < sub_frame_l/2-20) or \
                    (new_sources[min_distance_index]['ycentroid'] > sub_frame_l/2+20) or (new_sources[min_distance_index]['ycentroid'] < sub_frame_l/2-20):
                    
                    print('Large centroid deviation determined for target ', id, ' in image ',filename,'. Flagging.' )
                    save_centroiding_flags[id,i] = 1
                    if SEx_fwhm == max_SEx_fwhm:
                        print('Non-centered centroid determined for target ',id,', flagging.')
                        return ((cen_position[0],cen_position[1]),((sub_frame_l+1)/2,(sub_frame_l+1)/2),frame)
                    else: 
                        SEx_fwhm = SEx_fwhm + 1
                        continue
                
                else:
                    if id == 0:
                        #Estimate the seeing around the target
                        try:
                            (fwhm_x,fwhm_y) = seeing_estimator(frame,int(round(new_sources[min_distance_index]['xcentroid'])),int(round(new_sources[min_distance_index]['ycentroid'])))
                        except:
                            #If the gaussian fit fails, just assign the previous seeing measurement
                            print('Gaussian fit in seeing_estimator failed, assigning previous seeing values.')
                            fwhm_x = save_x_seeing[i-1]
                            fwhm_y = save_y_seeing[i-1]
                        save_x_seeing[i] = fwhm_x
                        save_y_seeing[i] = fwhm_y
                    #If a shift was applied correctly, should be within 2 pixels of the box center. 
                    if (abs(new_sources[min_distance_index]['xcentroid']-sub_frame_l/2) > 2.5) or ((abs(new_sources[min_distance_index]['ycentroid']-sub_frame_l/2) > 2.5) ):
                        #print('Non-centered centroid determined for id = ',str(id),', trying higher fwhm.')
                        if (SEx_fwhm == max_SEx_fwhm-1):
                            #pdb.set_trace()
                            save_centroiding_flags[id,i] = 1
                            return ((cen_position[0],cen_position[1]),((sub_frame_l+1)/2,(sub_frame_l+1)/2),frame)
                        else:
                            SEx_fwhm = SEx_fwhm + 1
                            continue
                    else:
                        if plot_subframes:
                            #Plot the measured centroid position. 
                            centroid1, = ax1.plot([new_sources[min_distance_index]['xcentroid']], [new_sources[min_distance_index]['ycentroid']], 'rx')
                            plt.pause(0.1)
                            centroid1.remove()
                            ax1.clear()
                        #if SEx_fwhm == max_SEx_fwhm-1:
                        #    pdb.set_trace()
                        return((new_targ_x,new_targ_y),(new_sources[min_distance_index]['xcentroid'],new_sources[min_distance_index]['ycentroid']), frame)

    def seeing_estimator(sub_frame,x_pos,y_pos):
        #Define Gaussian fit function
        def gaus(x,a,x0,sigma):
            return a*np.exp(-(x-x0)**2/(2*sigma**2))

        #First, do x dimension
        x_slice = sub_frame[int(sub_frame_l/2),:]-np.median(sub_frame[int(sub_frame_l/2),:])	
        x = np.arange(0,len(x_slice),1)        
        if i == 0:
            popt,pcov = curve_fit(gaus,x,x_slice,p0=[max(x_slice),int(sub_frame_l/2),8])
        else:
            popt,pcov = curve_fit(gaus,x,x_slice,p0=[max(x_slice),int(sub_frame_l/2),4])
        fwhm_x = popt[2]*2.2355*plate_scale #FWHM in arcsec
        if fwhm_x < 0:
            fwhm_x = abs(fwhm_x)

        #Now y
        y_slice = sub_frame[:,int(sub_frame_l/2)]-np.median(sub_frame[:,int(sub_frame_l/2)])
        if i == 0:
            popt,pcov = curve_fit(gaus,x,y_slice,p0=[max(y_slice),int(sub_frame_l/2),8])
        else:
            popt,pcov = curve_fit(gaus,x,y_slice,p0=[max(y_slice),int(sub_frame_l/2),4])

        fwhm_y = popt[2]*2.2355*plate_scale #FWHM in arcsec
        if fwhm_y < 0:
            fwhm_y = abs(fwhm_y)
        #print('X, Y Seeing (arcsec): ',np.round(fwhm_x,1),np.round(fwhm_y,1))

        if (fwhm_x > 5) or (fwhm_y > 5):
            print('Warning: Large seeing measured.')
            #pdb.set_trace()
            

        if (fwhm_x == 0) or (fwhm_y == 0):
            print('Warning: Zero seeing measured.')
            #pdb.set_trace()
        return(fwhm_x,fwhm_y)


    local_path = '/Users/tamburo/Documents/Data/PINES/'
    reduced_path = local_path+'Objects/'+target_name+'/domeflat_reduced/'
    reduced_files = natsort.natsorted(glob.glob(reduced_path+'*red.fits'))
    aper_phot_path = local_path+'Objects/'+target_name+'/aper_phot/'
    target_and_ref_path = aper_phot_path+'targ_and_refs/'
    image_shift_path = aper_phot_path+'image_shifts/'
    bpm = (1-pickle.load(open('/Users/tamburo/Documents/Data/PINES/Calibrations/Bad Pixel Masks/bpm.p','rb'))).astype('bool') 
    plate_scale = 0.579

    output_path = aper_phot_path+'centroids/'
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    #Check for existing output.
    if os.path.exists(output_path+target_name+'_centroids.p'):
        old_output = pickle.load(open(output_path+target_name+'_centroids.p','rb'))
        old_file_list = old_output['Old frame list']
        old_x_pos = old_output['X pos']
        old_y_pos = old_output['Y pos']
        old_x_subframe_pos = old_output['X subframe pos']
        old_y_subframe_pos = old_output['Y subframe pos']
        old_seeing_output = pickle.load(open(output_path+target_name+'_seeing.p','rb'))
        old_x_seeing = old_seeing_output['X seeing']
        old_y_seeing = old_seeing_output['Y seeing']
    else:
        old_file_list = []

    #Get list of files
    all_frame_list = np.array([])
    for i in range(np.size(reduced_files)):
        file_name = reduced_files[i].split('/')[-1]
        all_frame_list = np.append(all_frame_list, file_name)
    num_files = len(reduced_files)

    #Get initial positions of target and reference stars, output from identify_target_and_references.py.
    initial_positions = pickle.load(open(target_and_ref_path+target_name+'_target_and_ref_initial_locs.p','rb'))
    num_targets = len(initial_positions)

    #Get measured image shifts, output from deshifter.py.
    image_shifts = pickle.load(open(image_shift_path+target_name+'_image_shifts.p','rb'))
    image_x_shifts = image_shifts[0]
    image_y_shifts = image_shifts[1]

    #Set up figure for plotting the target/reference subframes. 
    if (plot_subframes):
        fig, ax1 = plt.subplots(1, 1, figsize=(6, 5))

    #Declare arrays for saving centroids as both aboslute positions and subframe positions.
    save_x_position = np.zeros((num_targets,num_files)) #The measured x position on the detector
    save_y_position = np.zeros((num_targets,num_files)) #The measured y position on the detector
    save_x_position_subframe = np.zeros((num_targets,num_files)) #The x position within the subframe
    save_y_position_subframe = np.zeros((num_targets,num_files)) #The y position within the subframe
    save_x_seeing = np.zeros(num_files) #Estimated seeing for the target in x direction
    save_y_seeing = np.zeros(num_files) #Estimated seeing for target in y direction
    save_centroiding_flags = np.zeros((num_targets,num_files)) #TODO: Flag bad centroiding. 
    #save_subframes = np.zeros((num_targets,num_files,sub_frame_l+1,sub_frame_l+1),dtype=np.float32)

    #Loop over all files
    first_hit = 0
    for i in range(num_files):
        filename = reduced_files[i].split('/')[-1]
        if filename not in old_file_list:
            print('Centroiding ',filename,', ',i+1,' of ',np.size(reduced_files))
            frame = fits.open(reduced_files[i])[0].data
            new_cen_positions = []

            if first_hit == 0: 
                prev_new_cen_positions = np.copy(initial_positions)
            
            first_hit += 1 

            #Declare some empty arrays into which we'll load subframes/masks/centroids of the target/references. 
            #sub_frame_stack = np.zeros((num_targets,sub_frame_l+1,sub_frame_l+1))
            #mask_stack = np.zeros((num_targets,sub_frame_l+1,sub_frame_l+1),dtype=bool)
            centroid_stack_x = np.zeros(num_targets)
            centroid_stack_y = np.zeros(num_targets)
            ap_flux_stack = np.zeros(num_targets)
            an_mean_stack = np.zeros(num_targets)
            avg3_stack = np.zeros(num_targets)
            med3_stack = np.zeros(num_targets)
            stddev3_stack = np.zeros(num_targets)

            #Use  measured shifts to accurately predict where the target/references will be in this frame. 
            for obj in range(len(prev_new_cen_positions)):
                prev_new_cen_positions[obj][0] = initial_positions[obj][0] + image_x_shifts[i]
                prev_new_cen_positions[obj][1] = initial_positions[obj][1] + image_y_shifts[i]
            
            #First, the target star.
            #Define a "sub_frame" around the target, to pass to both get_source and photometry. This speeds things up considerably 
            #	compared to passing the whole LMI image. This is done for the target, references, and cbv stars. 
            
            sub_frame = frame[int(prev_new_cen_positions[0][1])-sub_frame_l//2:int(prev_new_cen_positions[0][1]) + sub_frame_l//2+1, int(prev_new_cen_positions[0][0])-sub_frame_l//2:int(prev_new_cen_positions[0][0]) + sub_frame_l//2+1]
            
            #sub_frame = frame[int(prev_new_cen_positions[0][1])-sub_frame_l//2:int(prev_new_cen_positions[0][1]) + sub_frame_l//2+1, int(prev_new_cen_positions[0][0])-sub_frame_l//2:int(prev_new_cen_positions[0][0]) + sub_frame_l//2+1]
            bpm_frame = bpm[int(prev_new_cen_positions[0][1])-sub_frame_l//2:int(prev_new_cen_positions[0][1]) + sub_frame_l//2+1, \
                                int(prev_new_cen_positions[0][0])-sub_frame_l//2:int(prev_new_cen_positions[0][0]) + sub_frame_l//2+1]
            get_source_output = get_source(sub_frame,bpm_frame,prev_new_cen_positions[0],0)
            
            #sub_frame_stack[0,:,:] = get_source_output[2]
            #save_subframes[0,i,:,:] = np.uint16(get_source_output[2]) #Save the subframe as uint16 (to save memory). 
            new_cen_positions.append(get_source_output[0])
            save_x_position[0,i] = new_cen_positions[0][0]
            save_y_position[0,i] = new_cen_positions[0][1]
            centroid_stack_x[0] = get_source_output[1][0]
            centroid_stack_y[0] = get_source_output[1][1]
            save_x_position_subframe[0,i] = centroid_stack_x[0]
            save_y_position_subframe[0,i] = centroid_stack_y[0]

            #Then, all the references stars.
            for star in range(num_targets-1):
                sub_frame = frame[int(prev_new_cen_positions[star+1][1])-sub_frame_l//2:int(prev_new_cen_positions[star+1][1]) + sub_frame_l//2+1, \
                                    int(prev_new_cen_positions[star+1][0])-sub_frame_l//2:int(prev_new_cen_positions[star+1][0]) + sub_frame_l//2+1]
                bpm_frame = bpm[int(prev_new_cen_positions[star+1][1])-sub_frame_l//2:int(prev_new_cen_positions[star+1][1]) + sub_frame_l//2+1, \
                                    int(prev_new_cen_positions[star+1][0])-sub_frame_l//2:int(prev_new_cen_positions[star+1][0]) + sub_frame_l//2+1]
                #Correct the reference sub_frames for bad pixels.
                #sub_frame = ref_frame_corrector(sub_frame,bpm_frame)
                get_source_output = get_source(sub_frame,bpm_frame,prev_new_cen_positions[star+1],star+1)

                #sub_frame_stack[star+1,:,:] = get_source_output[2]
                #save_subframes[star+1,i,:,:] = np.uint16(get_source_output[2]) #Save the subframe as uint16 (to save memory). 
                try:
                    new_cen_positions.append(get_source_output[0])
                except:
                    pdb.set_trace()
                save_x_position[star+1,i] = new_cen_positions[star+1][0]
                save_y_position[star+1,i] = new_cen_positions[star+1][1]
                centroid_stack_x[star+1] = get_source_output[1][0]
                centroid_stack_y[star+1] = get_source_output[1][1]
                save_x_position_subframe[star+1,i] = centroid_stack_x[star+1]
                save_y_position_subframe[star+1,i] = centroid_stack_y[star+1]
            print('X, Y Seeing (arcsec): ',np.round(save_x_seeing[i],1),np.round(save_y_seeing[i],1))
            print('')
        #End loop over all files
        else:
            print(filename,' already has centroider output, skipping.')
            loc = np.where(np.array(old_file_list) == filename)[0][0]
            save_x_position[:,i] = old_x_pos[:,loc]
            save_y_position[:,i] = old_y_pos[:,loc]
            save_x_position_subframe[:,i] = old_x_subframe_pos[:,loc]
            save_y_position_subframe[:,i] = old_y_subframe_pos[:,loc]
            save_x_seeing[i] = old_x_seeing[loc]
            save_y_seeing[i] = old_y_seeing[loc]

    #Print out the number of files for each target with bad centroiding, might help you through out certain references. 
    print('')
    print('% of files with bad centroiding')
    print('Target: ', np.round(100*len(np.where(save_centroiding_flags[0,:] == 1)[0])/num_files,1))
    for ref in range(num_targets-1):
        print('Ref. ',str(ref),': ',np.round(100*len(np.where(save_centroiding_flags[ref+1,:] == 1)[0])/num_files,1))

    centroid_output_dict = {'X pos':save_x_position,'Y pos':save_y_position,'X subframe pos':save_x_position_subframe,
                            'Y subframe pos':save_y_position_subframe,'Old frame list':all_frame_list}

    #Save positions and seeing as pickle files.  
    pickle.dump(centroid_output_dict,open(output_path+target_name+'_centroids.p','wb'))
    
    seeing_output_dict = {'X seeing':save_x_seeing,'Y seeing':save_y_seeing}
    
    pickle.dump(seeing_output_dict,open(output_path+target_name+'_seeing.p','wb'))
    
    print('Done centroiding!')
    print('')
