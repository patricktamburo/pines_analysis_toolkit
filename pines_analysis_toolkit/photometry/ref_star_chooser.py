from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.photometry.detect_sources import detect_sources
from pines_analysis_toolkit.photometry.target_finder import target_finder
from pines_analysis_toolkit.utils.quick_plot import quick_plot
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from photutils import CircularAperture, aperture_photometry, CircularAnnulus
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from scipy.stats import sigmaclip
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader
import time
from natsort import natsorted

def ref_star_chooser(short_name, source_detect_image, guess_position=(700.,382.), radius_check=6., non_linear_limit=3300., 
                    dimness_tolerance=0.5, brightness_tolerance=10., closeness_tolerance=12., distance_from_target=900., edge_tolerance=50., exclude_lower_left=False, restore=False,
                    source_detect_plot=False, force_output_path=''):
    """Chooses suitable reference stars for a target in a specified source_detect_image.

    :param short_name: the short name for the target
    :type short_name: str
    :param source_detect_image: name of the reduced image in which you want to find reference stars
    :type source_detect_image: str
    :param guess_position: guess position for the target in the source_detect_image, defaults to (700.,382.)
    :type guess_position: tuple, optional
    :param radius_check: the radius in pixels used to perform photometry to compare target and reference brightnesses, defaults to 6.
    :type radius_check: float, optional
    :param non_linear_limit: ADU value above which references are considered to be in the non-linear limit of the detecto, defaults to 3300.
    :type non_linear_limit: float, optional
    :param dimness_tolerance: minimum multiple of the target's measured brightness that a reference is allowed to have, defaults to 0.5
    :type dimness_tolerance: float, optional
    :param brightness_tolerance: maximum multiple of the target's measured brightness that a reference is allowed to have, defaults to 10.
    :type brightness_tolerance: float, optional
    :param closeness_tolerance: the closest distance in pixels that a reference star can be to another detected source and still be considered as a reference, defaults to 12.
    :type closeness_tolerance: float, optional
    :param distance_from_target: the furthest distance in pixels that a reference star can be from the target and still be considered, defaults to 900.
    :type distance_from_target: float, optional
    :param edge_tolerance: the closest distance in pixels that a reference can be to the edge of the detector and still be considered, defaults to 50.
    :type edge_tolerance: float, optional
    :param exclude_lower_left: whether or not to exclude reference stars from the lower left quadrant (due to occasional Mimir 'bars' issue), defaults to False
    :type exclude_lower_left: bool, optional
    :param restore: whether or not to restore references from previous output, defaults to False
    :type restore: bool, optional
    :param source_detect_plot: whenther or not to plot all detected sources, defaults to False
    :type source_detect_plot: bool, optional
    :param force_output_path: top-level 'force output path', used if you want to use a folder other than ~/Documents/PINES_analysis_toolkit/, defaults to ''
    :type force_output_path: str, optional
    :raises ValueError: If the measured seeing in your selected source_detect_image is Nan
    :raises ValueError: If the measured x/y shift in your selected_source_detect_image is > 20 pixels
    :return: Saves a plot of the target/detected reference stars to the object's 'sources' directory. Saves .csv of target/reference pixel positions in the object's 'sources' directory.
    :rtype: plot/csv
    """

    plt.ion()
    #Get your local PINES directory
    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()

    if restore:
        output_df = pd.read_csv(pines_path/('Objects/'+short_name+'/sources/target_and_references_source_detection.csv'))
        print('')
        print('Restoring ref_star_chooser output from {}.'.format(pines_path/('Objects/'+short_name+'/sources/target_and_references_source_detection.csv')))
        print('')
        return output_df

    #Find reduced files in the directory.
    data_dir = pines_path/('Objects/'+short_name+'/reduced/')
    reduced_files = np.array(natsorted([x for x in data_dir.glob('*.fits')]))

    #If the directory doesn't exist, or there are no reduced files, return nothing. 
    if (not data_dir.exists()) or (len(reduced_files) == 0):
        print('ERROR: No reduced images exist for {}.'.format(short_name))
        return

    #Set path to source directory. 
    source_dir = pines_path/('Objects/'+short_name+'/sources/')

    source_detect_image_ind = np.where([i.name == source_detect_image for i in reduced_files])[0][0]
    source_detect_image_path = reduced_files[source_detect_image_ind]
    source_frame = source_detect_image_path.name.split('_')[0]+'.fits'
    log_name = source_frame.split('.')[0]+'_log.txt'
    log = pines_log_reader(pines_path/('Logs/'+log_name))
    source_frame_ind = np.where(log['Filename'] == source_frame)[0][0]
    extra_x_shift = log['X shift'][source_frame_ind]
    extra_y_shift = log['Y shift'][source_frame_ind]
    source_detect_seeing = float(log['X seeing'][source_frame_ind])

    #Make sure the seeing value isn't a NaN. 
    if np.isnan(source_detect_seeing):
        raise ValueError('Seeing value = NaN in log. Try a different source_detect_image_ind in call to ref_star_chooser.')
    
    if (float(extra_x_shift) > 20) or (float(extra_y_shift) > 20):
        raise ValueError('Measured x or y shift > 20 pixels. Try a different source_detect_image_ind in call to ref_star_chooser.')
    
    extra_shift_path = (source_dir/'extra_shifts.txt')
    with open(extra_shift_path, 'w') as file:
        file.write(str(extra_x_shift)+' '+str(extra_y_shift))

    #source_detect_seeing = 3.5

    #Detect sources in the image. 
    sources = detect_sources(source_detect_image_path, source_detect_seeing, edge_tolerance, plot=source_detect_plot, thresh=3.0)

    #Identify the target in the image using guess_position. 
    target_id = target_finder(sources, guess_position)

    #Now determine ids of suitable reference stars.
    image = fits.open(source_detect_image_path)[0].data

    #Interpolate nans if any in image. 
    kernel = Gaussian2DKernel(x_stddev=0.5)
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
            brightness_flag = ((aperture_photometry(image, ap)['aperture_sum'] - bg_estimate * ap.area) < brightness_tolerance*targ_flux_estimate)[0]
            closeness_flag = (len(np.where(dists[np.where(dists != 0)] < closeness_tolerance)[0]) == 0)
            proximity_flag = (np.sqrt((sources['xcenter'][target_id]-sources['xcenter'][i])**2+(sources['ycenter'][target_id]-sources['ycenter'][i])**2) < distance_from_target)
            if exclude_lower_left:
                if (sources['xcenter'][i] < 512) & (sources['ycenter'][i] < 512):
                    lower_left_flag = False
                else:
                    lower_left_flag = True
            else:
                lower_left_flag = True

            if (non_linear_flag) & (dimness_flag) & (brightness_flag) & (closeness_flag) & (proximity_flag) & (lower_left_flag):
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
    im = ax.imshow(image, vmin=stats[1], vmax=stats[1]+7*stats[2], origin='lower', cmap='Greys')
    cax = fig.add_axes([0.87, 0.15, 0.035, 0.7])
    fig.colorbar(im, cax=cax, orientation='vertical', label='ADU')
    ax.set_title(source_detect_image_path.name)
    ax.plot(output_df['Source Detect X'][0], output_df['Source Detect Y'][0], marker='o', color='m', mew=2, ms=12, ls='', mfc='None', label='Target')
    #Plot selected references
    for i in range(1,len(output_df)):
        ax.plot(output_df['Source Detect X'][i], output_df['Source Detect Y'][i], marker='o', color='b', mew=2, ms=12, ls='', mfc='None', label='Reference')
        ax.text(output_df['Source Detect X'][i]+7, output_df['Source Detect Y'][i]+5, str(i), fontsize=14, color='r')
    #Plot detected sources
    for i in range(len(sources)):
        ax.plot(sources['xcenter'], sources['ycenter'], 'rx', label='Detected sources', alpha=0.6, ms=3)

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
            ax.text(output_df['Source Detect X'][i]+7, output_df['Source Detect Y'][i]+5, str(i), color='r', fontsize=14)
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
        print('Saving target and references image to {}.'.format(source_dir/('target_and_refs.png')))
        plt.savefig(source_dir/('target_and_refs.png'))
        plt.close('all')
        return output_df
    elif ans =='n':
        raise ValueError('Try changing arguments of ref_star_chooser or detect_sources and try again.')
    else:
        return