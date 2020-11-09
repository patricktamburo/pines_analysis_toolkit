from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.photometry.detect_sources import detect_sources
from pines_analysis_toolkit.photometry.target_finder import target_finder
from pines_analysis_toolkit.utils.quick_plot import quick_plot
import pdb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from photutils import CircularAperture, aperture_photometry, CircularAnnulus
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from scipy.stats import sigmaclip
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader
import time

'''Authors:
		Patrick Tamburo, Boston University, June 2020
	Purpose:
        Chooses suitable reference stars for a target.
	Inputs:
        target (str): The target's full 2MASS name
        source_detect_image (str, optional): The name of the reduced file to use for source detection. By default, the program will use the master image for the target. 
        seeing_fwhm (float, optional): Seeing in arcsec for the source detection.
        radius_check (float, optional): The radius (in pixels) of an aperture that will be placed around potential reference stars to see if there are non-linear pixels.
        non_linear_limit (float, optional): If a potential reference star has a pixel within an aperture of radius radius_check with a value greater than non_linear_limit, it will not be used.
            NOTE: Mimir data becomes non-linear around 4000 ADU.
        dimness_tolerance (float, optional): The minimum fraction of the target's brightness that a reference star may have. 
        closeness_tolerance (float, optional): The closest a potential reference star can be to another source before it's discarded (in pixels). 
        distance_from_target (float, optional): The furthest a potential reference star can be from the target before it's discarded (in pixels).
        exclude_lower_left (bool, optional): Whether or not to automatically discard potential references in the lower left quadrant, which is necessary if Mimir's bars problem exists in the data. 
        restore (bool, optional): Whether or not to restore ref_star_chooser output that already exists. 

    Outputs:
		ref stars...
	TODO:
'''

def ref_star_chooser(target, source_detect_image='', guess_position=(705.,386.), seeing_fwhm=2.5, radius_check=6., non_linear_limit=3300., dimness_tolerance=0.5, closeness_tolerance=12., distance_from_target=700., exclude_lower_left=False, restore=False):
    plt.ion() 
    
    #Get the target's 'short name'
    short_name = short_name_creator(target)

    #Get your local PINES directory
    pines_path = pines_dir_check()

    if restore:
        output_df = pd.read_csv(pines_path/('Objects/'+short_name+'/sources/target_and_references_source_detection.csv'))
        print('')
        print('Restoring ref_star_chooser output from {}.'.format(pines_path/('Objects/'+short_name+'/sources/target_and_references_source_detection.csv')))
        print('')
        return output_df

    #Find reduced files in the directory.
    data_dir = pines_path/('Objects/'+short_name+'/reduced/')
    reduced_files = [x for x in data_dir.glob('*.fits')]

    #If the directory doesn't exist, or there are no reduced files, return nothing. 
    if (not data_dir.exists()) or (len(reduced_files) == 0):
        print('ERROR: No reduced images exist for {}.'.format(target))
        return

    #Set path to source directory. 
    source_dir = pines_path/('Objects/'+short_name+'/sources/')

    #By default, use the target's master image for the source detection. 
    if source_detect_image == '':
        source_detect_image_path = pines_path/('Master Images/'+target.replace(' ','')+'_master.fits')
        #If the master image is used for source detection, centroids will be where we expect them to be based on the measured x/y shifts in the log,
        #and we don't need to apply any "extra" x/y shifts. 
        extra_x_shift = '0.0'
        extra_y_shift = '0.0'
    #If a different reduced file is specified, use that. 
    else:
        source_detect_image_path = data_dir/(source_detect_image)
        source_frame = source_detect_image.split('_')[0]+'.fits'
        log_name = source_frame.split('.')[0]+'_log.txt'
        log = pines_log_reader(pines_path/('Logs/'+log_name))
        source_frame_ind = np.where(log['Filename'] == source_frame)[0][0]
        #If the master image is *not* used for source detection, the image will shifted with respect to the master image, and extra x/y shifts must be
        #applied in order for the x/y shifts in the log to accurately predict where the sources will be. 
        extra_x_shift = log['X shift'][source_frame_ind]
        extra_y_shift = log['Y shift'][source_frame_ind]
    
    extra_shift_path = (source_dir/'extra_shifts.txt')
    with open(extra_shift_path, 'w') as file:
        file.write(str(extra_x_shift)+' '+str(extra_y_shift))

    #Detect sources in the image. 
    sources = detect_sources(source_detect_image_path, source_dir, seeing_fwhm, thresh=2., plot=True)

    target_id = target_finder(sources, guess_position)

    #Now determine ids of suitable reference stars.
    image = fits.open(source_detect_image_path)[0].data

    #Interpolate nans if any in image. 
    kernel = Gaussian2DKernel(x_stddev=1)
    image = interpolate_replace_nans(image, kernel)

    suitable_ref_ids = []
    bad_ref_ids = [target_id]
    #Get a value for the brightness of the target with a radius_check aperture. We don't want reference stars dimmer than dimness_tolerance * targ_flux_estimates. 
    target_ap = CircularAperture((sources['xcenter'][target_id], sources['ycenter'][target_id]), r=radius_check)
    target_an = CircularAnnulus((sources['xcenter'][target_id], sources['ycenter'][target_id]), r_in=12, r_out=30)
    mask = target_an.to_mask(method='center')
    an_data = mask.multiply(image)
    an_data_1d = an_data[an_data != 0]
    vals, lower, upper = sigmaclip(an_data_1d)
    bg_estimate = np.median(vals)
    targ_flux_estimate = aperture_photometry(image, target_ap)['aperture_sum'] - bg_estimate * target_ap.area
    
    for i in range(len(sources)):
        if i != target_id:
            potential_ref_loc = (sources['xcenter'][i], sources['ycenter'][i])
            #Identify brightnesses of pixels surrounding the potential reference within a 10 pixel radius. 
            ap = CircularAperture(potential_ref_loc, r=radius_check)
            an = CircularAnnulus(potential_ref_loc, r_in=12, r_out=30)
            mask = an.to_mask(method='center')
            an_data = mask.multiply(image)
            an_data_1d = an_data[an_data != 0]
            vals, lower, upper = sigmaclip(an_data_1d)
            bg_estimate = np.median(vals)

            pixels = ap.to_mask()
            sub_image = pixels.cutout(image)

            dists = []
            for j in range(len(sources)):
                dists.append(np.sqrt((sources['xcenter'][i]-sources['xcenter'][j])**2+(sources['ycenter'][i]-sources['ycenter'][j])**2))
            dists = np.array(dists)
            
            non_linear_flag = len(np.where(sub_image > non_linear_limit)[0]) == 0 #Check if potential ref is near non-linear regime.
            dimness_flag = ((aperture_photometry(image, ap)['aperture_sum'] - bg_estimate * ap.area) > dimness_tolerance*targ_flux_estimate)[0] #Check if potential ref is bright enough. 
            closeness_flag = (len(np.where(dists[np.where(dists != 0)] < closeness_tolerance)[0]) == 0)
            proximity_flag = (np.sqrt((sources['xcenter'][target_id]-sources['xcenter'][i])**2+(sources['ycenter'][target_id]-sources['ycenter'][i])**2) < distance_from_target)
            if exclude_lower_left:
                if (sources['xcenter'][i] < 512) & (sources['ycenter'][i] < 512):
                    lower_left_flag = False
                else:
                    lower_left_flag = True
            else:
                lower_left_flag = True

            if (non_linear_flag) & (dimness_flag) & (closeness_flag) & (proximity_flag) & (lower_left_flag):
                suitable_ref_ids.append(i)
            else:
                bad_ref_ids.append(i)

    #ax.plot(sources['xcenter'][suitable_ref_ids], sources['ycenter'][suitable_ref_ids],'yx')
    print('Found {} suitable reference stars.'.format(len(suitable_ref_ids)))
    print('')
    
    if len(suitable_ref_ids) < 2:
        print('Not enough reference stars found. Try loosening reference star criteria.')
        return
    
    target = sources.iloc[[target_id]].reset_index(drop=True).drop(columns=['id'])
    suitable_refs = sources.drop(bad_ref_ids).reset_index(drop=True).drop(columns=['id'])

    output_df = pd.DataFrame(columns = ['Name','Source Detect X','Source Detect Y'])
    output_df = output_df.append({'Name': short_name, 'Source Detect X': sources['xcenter'][target_id], 'Source Detect Y': sources['ycenter'][target_id]}, ignore_index=True)
    for i in range(len(suitable_refs)):
        output_df = output_df.append({'Name': 'Reference '+str(i+1), 'Source Detect X': sources['xcenter'][suitable_ref_ids[i]], 'Source Detect Y': sources['ycenter'][suitable_ref_ids[i]]}, ignore_index=True)
    
    
    stats = sigma_clipped_stats(image)
    fig, ax = plt.subplots(figsize=(10,9))
    ax.set_aspect('equal')
    im = ax.imshow(image, origin='lower', vmin=stats[1], vmax=stats[1]+7*stats[2], cmap='Greys')
    cax = fig.add_axes([0.87, 0.15, 0.035, 0.7])
    fig.colorbar(im, cax=cax, orientation='vertical', label='ADU')
    ax.set_title(source_detect_image_path.name)
    ax.plot(output_df['Source Detect X'][0], output_df['Source Detect Y'][0], marker='o', color='m', mew=2, ms=12, ls='', mfc='None', label='Target')
    for i in range(1,len(output_df)):
        ax.plot(output_df['Source Detect X'][i], output_df['Source Detect Y'][i], marker='o', color='b', mew=2, ms=12, ls='', mfc='None', label='Reference')
        ax.text(output_df['Source Detect X'][i]+5, output_df['Source Detect Y'][i]+5, str(i), fontsize=14, color='r')
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.1, 1.1))
    plt.show()

    #Sometimes, the source detection returns things that are clearly not reference stars. 
    #Allow the user to remove them here. 
    ans = input('Enter IDs of references to remove, separated by commas (e.g.: 1,4,8,14).\nIf none to remove, hit enter:  ')
    if ans != '':
        ans = ans.split(',')
        ref_ids_to_remove = []
        for i in range(len(ans)):
            ref_ids_to_remove.append(int(ans[i]))

        ref_ids_to_remove = np.sort(ref_ids_to_remove)[::-1] #Reverse sort the list. 
        #Drop those references from output_df
        for i in ref_ids_to_remove:
            output_df = output_df.drop(index=i)
        output_df = output_df.reset_index(drop=True)
        #Rename the remaining references. 
        for i in range(1,len(output_df)):
            name = output_df['Name'][i]
            if int(name.split(' ')[1]) != i:
                name = 'Reference '+str(i)
                output_df['Name'][i] = name
        #Replot everything. 
        ax.cla()
        im = ax.imshow(image, origin='lower', vmin=stats[1], vmax=stats[1]+7*stats[2], cmap='Greys')
        cax = fig.add_axes([0.87, 0.15, 0.035, 0.7])
        fig.colorbar(im, cax=cax, orientation='vertical', label='ADU')
        ax.set_title(source_detect_image_path.name)
        ax.plot(output_df['Source Detect X'][0], output_df['Source Detect Y'][0], marker='o', color='m', mew=2, ms=12, ls='', mfc='None', label='Target')
        for i in range(1,len(output_df)):
            ax.plot(output_df['Source Detect X'][i], output_df['Source Detect Y'][i], marker='o', color='b', mew=2, ms=12, ls='', mfc='None', label='Reference')
            ax.text(output_df['Source Detect X'][i]+5, output_df['Source Detect Y'][i]+5, str(i), fontsize=14)
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.1, 1.1))
    plt.show()

    ans = input('Happy with reference star selection? y/n: ')
    if ans == 'y':
        print('')
        print('Saving target/reference info to {}.'.format(source_dir/'target_and_references_source_detection.csv'))
        output_df.to_csv(source_dir/'target_and_references_source_detection.csv',index=False)
        print('')
        print('Saving target and references image to {}.'.format(source_dir/(source_detect_image_path.name.split('.fits')[0]+'_target_and_refs.png')))
        plt.savefig(source_dir/(source_detect_image_path.name.split('.fits')[0]+'_target_and_refs.png'))
        plt.close('all')
        return output_df
    elif ans =='n':
        raise ValueError('Try changing arguments of ref_star_chooser or detect_sources and try again.')
    else:
        return