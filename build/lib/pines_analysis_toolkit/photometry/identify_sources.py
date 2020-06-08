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
import glob
import natsort

def identify_sources(target_name,target_location_guess,source_frame_id=0,SEx_fwhm=8.,SEx_sigma=3.,
                                   edge_cut=50.,star_star_dist=20.,snr_cut=0.25,non_linear_cut=4000.,distance_cut=700.):
    
    """ 
        AUTHOR: 
            Patrick Tamburo, BU, Dec. 2019
        PURPOSE: 
            Identifies your target in a Mimir image, and finds suitable reference stars. 
        INPUTS:
            target_name (str): The name of your target, which must match the data analysis folder name. 
            local_path (str): The analysis directory for this target.
            target_location_guess (tuple of ints): Your guess for the target's position in the source_detect_frame. 
            SEx_fwhm (float): The fwhm used in DAOStarFinder to detect stars. By default, set to 4. 
            SEx_sigma (float): The sigma above background used in DAOStarFinder to detect stars. By default, set to 3. 
            edge_cut (float): The minimum distance a star must have from the edge of the detector to be considered as a reference.
                By default, set to 50. 
            star_star_dist (float): The minimum distance a star must have from other identified sources in order to be considered
                as a reference. By default, set to 20. 
            snr_cut (flaot): The fraction of the target SNR that a star must have to be considered a reference. Values between 0 and 1. 
                By default, set to 0.5. 
            non_linear_cut (float): The peak counts below which a star is considered linear, and used as a potential reference. 
                TODO: figure out what this value should be for Mimir. I think the default value of 4000 is too high. 
            distance_cut (float): The maximum distance allowed between the target and potential reference stars. By default, set to 700. 
        OUTPUTS:
            Writes four files to local_path/aper_phot/targ_and_refs:
                1) target_name_identify_target_and_references_log.txt: A text log of your inputs and discovered sources/references. 
                2) target_name_sources_initial_locs.p: A pickle file containing the locations of all sources in the image. 
                    These are used to create a 'master_synthetic_image' in deshifter.py, which we use to calculate shifts between images.
                3) target_name_target_and_ref_initial_locs.p: A pickle file containing the positions of your target and reference stars in 
                    the source_detect_frame. These are used in centroider.py to determine the positions of all stars of interest in every image. 
                4) target_name_target_and_refs.png: An image of the sources detected in the image, your target, and suitable reference stars.
    """

    local_path = '/Users/tamburo/Documents/Data/PINES/'
    reduced_path = local_path+'Objects/'+target_name+'/domeflat_reduced/'
    aper_phot_path = local_path+'Objects/'+target_name+'/aper_phot/'
    bpm = (1-pickle.load(open('/Users/tamburo/Documents/Data/PINES/Calibrations/Bad Pixel Masks/bpm.p','rb'))).astype('bool') 
    output_path = aper_phot_path+'targ_and_refs/'

    if not os.path.exists(output_path):
        os.mkdir(output_path)
    
    reduced_files = natsort.natsorted(glob.glob(reduced_path+'*red.fits'))

    #Now open up the source_detect_frame, the first in the directory. 
    source_frame = fits.open(reduced_files[source_frame_id])[0].data
    frame_size = np.shape(source_frame)

    #Get sigma_clipped_stats of frame. 
    avg, med, stddev = sigma_clipped_stats(source_frame, sigma=3.0,maxiters=3)    
    print('Average, median, stddev.')
    print(avg, med, stddev)
    print('')

    #Find stars in the image. 
    #daofind = DAOStarFinder(fwhm=SEx_fwhm, threshold=SEx_sigma*stddev, ratio=0.9, exclude_border=True, sky=np.nanmedian(source_frame-med))
    daofind = DAOStarFinder(fwhm=SEx_fwhm, threshold=SEx_sigma*stddev, ratio=0.9, exclude_border=True, sky=med)

    daosources = daofind(source_frame, mask=bpm)

    #Grab the source parmeters as arrays, I find that easier than working with the table.
    x_centroids = daosources['xcentroid']
    y_centroids = daosources['ycentroid']
    sharpness = daosources['sharpness']
    fluxes = daosources['flux']
    peaks = daosources['peak']

    #Cut sources that are found within 20 pix of the edges.
    use_pos = np.where((x_centroids > 20) & (x_centroids < 1004) & (y_centroids  > 20) & (y_centroids  < 1004))[0]
    x_centroids = x_centroids[use_pos]
    y_centroids = y_centroids[use_pos]
    sharpness = sharpness[use_pos]
    fluxes = fluxes[use_pos]
    peaks = peaks[use_pos]

    #Also cut on sharpness, this seems to eliminate a lot of false detections.
    use_sharp = np.where(sharpness > 0.4)[0]
    x_centroids  = x_centroids [use_sharp]
    y_centroids  = y_centroids [use_sharp]
    sharpness = sharpness[use_sharp]
    fluxes = fluxes[use_sharp]
    peaks = peaks[use_sharp]

    #Finally, cut targets whose y centroids are near y = 512. These are usually bad.
    use_512 = np.where(np.logical_or((y_centroids < 510),(y_centroids > 514)))[0]
    x_centroids  = x_centroids [use_512]
    y_centroids  = y_centroids [use_512]
    sharpness = sharpness[use_512]
    fluxes = fluxes[use_512]
    peaks = peaks[use_512]

    num_sources = len(x_centroids)
    print(num_sources,' sources found.')
    print()

    #Plot the sources found. 
    norm = simple_norm(source_frame, 'sqrt')
    plt.ion()
    fig = plt.figure(figsize=(10,9))
    ax = fig.add_subplot(1,1,1)
    im = ax.imshow(source_frame,origin='lower',vmin=med,vmax=med+7*stddev,cmap='Greys_r',norm=norm)
    fig.colorbar(im,label='Counts')
    ax.plot(x_centroids,y_centroids,'bo',mfc='None',ms=10,mew=2,alpha=0.8,label='Sources')
    for i in range (len(x_centroids)):
        ax.text(x_centroids[i]+5,y_centroids[i]+5,str(i),fontsize=14,color='b')
    plt.tight_layout()
    plt.show()

    print('Inspect the sources. If any need to be removed, note their ids. When you have identified them (or if there arent any), press c to continue.')
    print('')
    pdb.set_trace()

    print('')
    #If any false detections persisted through the cuts, identify them by their id and remove them. This is needed to accurately determine shifts between this frame
    #   and subsequent frames in your data, which we do by cross-correlating synthetic images of the fields based on source detections. 
    ids = input('Enter ids of sources to be removed separated by commas (i.e., 4,18,22). If none to remove, hit enter:  ')
    if ids != '':
        ids_to_eliminate = [int(i) for i in ids.split(',')]
        x_centroids  = np.delete(x_centroids,ids_to_eliminate)
        y_centroids  = np.delete(y_centroids,ids_to_eliminate)
        sharpness = np.delete(sharpness,ids_to_eliminate)
        fluxes = np.delete(fluxes,ids_to_eliminate)
        peaks = np.delete(peaks,ids_to_eliminate)
        plt.cla()
        im = ax.imshow(source_frame,origin='lower',vmin=med,vmax=med+7*stddev,cmap='Greys_r',norm=norm)
        ax.plot(x_centroids,y_centroids,'bo',mfc='None',ms=10,mew=2,alpha=0.8,label='Sources')
        for i in range (len(x_centroids)):
            ax.text(x_centroids[i]+5,y_centroids[i]+5,str(i),fontsize=14,color='b')
        plt.tight_layout()
        plt.show()
        num_sources = len(x_centroids)
        print(num_sources,' sources remain after cut.')
        print()

    #Now, find the target. 
    source_distances_from_guess = np.sqrt((x_centroids-target_location_guess[0])**2+(y_centroids-target_location_guess[1])**2)
    target_id = np.where(source_distances_from_guess == np.min(source_distances_from_guess))[0][0]
    print('The closest source found to your guess was id = ',target_id,', at (',np.round(x_centroids[target_id],1),',',np.round(y_centroids[target_id],1),').')
    print()

    #Add it to the plot.
    ax.plot(x_centroids[target_id],y_centroids[target_id],'ro',mfc='None',ms=10,mew=2,alpha=0.8,label='Target')
    ax.text(x_centroids[target_id]+5,y_centroids[target_id]+5,str(target_id),fontsize=14,color='r')

    #Now, find suitable reference stars. 
    print('Finding suitable reference stars.')

    #Zeroth criteria: the target can't be a reference.
    zeroth_cut = np.ones(num_sources,dtype=bool)
    zeroth_cut[target_id] = False

    #First criteria: Edge. The target centroid must be at least edge_cut pixels away from the edge of the detector.
    first_cut_x = np.ones(num_sources,dtype=bool)
    first_cut_y = np.ones(num_sources,dtype=bool)
    first_cut_x[np.where((x_centroids>(frame_size[1]-(edge_cut))) | (x_centroids<(edge_cut)))] = False
    first_cut_y[np.where((y_centroids>(frame_size[0]-(edge_cut))) | (y_centroids<(edge_cut)))] = False

    #Second criteria: Distance. Target centroid cannot be within star_star_dist of any other target. This may be loosened in the future for sparse fields.
    # Put all centroids in a Nx2 coord array, then run cdist to get distances between.
    source_coords = np.transpose(np.array([x_centroids,y_centroids]))
    distance_array = distance.cdist(source_coords,source_coords,'euclidean')
    second_cut = np.ones(num_sources,dtype=bool)
    #for i in range(size(source_coords)/2):
    for i in range(int(np.size(source_coords)/2)): 
        if np.sort(distance_array[i])[1] < star_star_dist:
            second_cut[i] = False
    num_eliminated_second_cut = len(np.where(second_cut == 0)[0])

    #Third: Saturation. Immediately exclude any target with peak counts greater than Nk.
    #TODO: Not sure if this appropriately identifies saturated sources for Mimir.
    third_cut = np.ones(num_sources,dtype=bool)
    third_cut[np.where(peaks>non_linear_cut)] = False
    num_eliminated_third_cut = len(np.where(third_cut == 0)[0])

    #Fourth: Brightness. Stars must have a S/N ratio of at least N% of the target. S/N
    # will be judged as the peak brightness divided by the median background
    fourth_cut = np.ones(num_sources,dtype=bool)
    snr_target = peaks[target_id]/med
    fourth_cut[np.where((np.array(peaks)/med)<(snr_cut*snr_target))] = False
    num_eliminated_fourth_cut = len(np.where(fourth_cut == 0)[0])

    #Fifth: Distance from target. Our fields are so crowded that we can take only a subset of
    #	references that are nearby. This will help minimize flat-fielding differences across
    #	the chip as well. 
    fifth_cut = np.ones(num_sources,dtype=bool)	
    fifth_cut[np.where(np.sqrt((x_centroids[target_id]-x_centroids)**2+(y_centroids[target_id]-y_centroids)**2) > distance_cut)] = False

    #TODO: Add a cut for references near bad pixels? 

    #Now set the ref_ids based on all the cuts.
    sources_id = np.arange(0,num_sources,1)
    ref_ids = sources_id[(zeroth_cut*first_cut_x*first_cut_y*second_cut*third_cut*fourth_cut*fifth_cut)]

    print('Found ',len(ref_ids),' suitable reference stars.')
    if len(ref_ids)  <= 2:
        print('Not enough reference stars. Try loosening reference star cut criteria.')
        pdb.set_trace()
    if len(ref_ids) > 10:
        #Take the 10 brighest
        print('Taking the 10 brightest suitable reference stars.')
        bright_locs = np.argsort(peaks[ref_ids])[::-1][0:10]
        ref_ids = ref_ids[bright_locs]
    ax.plot(x_centroids[ref_ids],y_centroids[ref_ids],'yo',mfc='None',ms=10,mew=2,label='References')
    for ref_id in ref_ids:
        ax.text(x_centroids[ref_id]+5,y_centroids[ref_id]+5,str(ref_id),fontsize=14,color='y')
    ax.legend(loc='upper left')
    plt.pause(0.5)
    print('')

    #Save the target and reference locations in a pickle file. By our convention, the target is the first object in these arrays. 
    save_centroids = []
    for i in range(len(ref_ids)+1):
        if i == 0:
            save_centroids.append((x_centroids[target_id],y_centroids[target_id]))
        else:
            save_centroids.append((x_centroids[ref_ids[i-1]],y_centroids[ref_ids[i-1]]))
    save_centroids = np.array(save_centroids)

    #Also save ALL of the source centroids, to make a synthetic image against which we'll cross-correlate to determine shifts of other images of this field. s
    save_source_centroids = []
    for i in range(num_sources):
        save_source_centroids.append((x_centroids[i],y_centroids[i]))
    save_source_centroids = np.array(save_source_centroids)

    pdb.set_trace()
    write_out = False
    output_filename = output_path+target_name+'_target_and_ref_initial_locs.p'
    if os.path.exists(output_filename):
        yn = input('Output file already exists, do you want to overwrite it? y/n:   ')
        if yn.lower() == 'y':
            write_out = True
        else:
            print('Output not saved.')
    else:
        write_out = True

    if write_out == True:
        print('Saving figure and output.')
        print('')
        plt.savefig(output_path+target_name+'_target_and_refs.png')
        pickle.dump(save_centroids,open(output_path+target_name+'_target_and_ref_initial_locs.p','wb'))
        pickle.dump(save_source_centroids,open(output_path+target_name+'_sources_initial_locs.p','wb'))

        #Write to an output_log file.
        output_log_filename = output_path+target_name+'_identify_target_and_references_log.txt'
        output_log = open(output_log_filename,'w')
        output_log.write('SOURCE DETECTION INPUTS\n')
        output_log.write('Target = '+target_name+'\n')
        output_log.write('Source detect frame = '+reduced_files[source_frame_id]+'\n')
        output_log.write('Target location guess = '+str(target_location_guess)+'\n')
        output_log.write('Source extraction fwhm = '+str(SEx_fwhm)+'\n')
        output_log.write('Source extraction sigma = '+str(SEx_sigma)+'\n')
        output_log.write('\n')
        output_log.write('REFERENCE STAR CUTS\n')
        output_log.write('Edge cut = '+str(edge_cut)+'\n')
        output_log.write('Star-star-distance cut = '+str(star_star_dist)+'\n')
        output_log.write('SNR cut = '+str(snr_cut)+'\n')
        output_log.write('Non-linear cut = '+str(non_linear_cut)+'\n')
        output_log.write('Distance cut = '+str(distance_cut)+'\n')
        output_log.write('\n')
        output_log.write('TARGET AND REFERENCE LOCATIONS\n')
        output_log.write('Number of sources found = '+str(num_sources)+'\n')
        output_log.write('Number of references found = '+str(len(ref_ids))+'\n')
        for i in range(len(save_centroids)):
            if i == 0:
                output_log.write('Target = '+str(save_centroids[i])+'\n')
            else:
                output_log.write('Ref. '+str(i)+' = '+str(save_centroids[i])+'\n')
        output_log.close()

    print('Done!')