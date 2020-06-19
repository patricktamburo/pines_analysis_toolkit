from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.photometry.detect_sources import detect_sources
from pines_analysis_toolkit.photometry.target_finder import target_finder
import pdb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from photutils import CircularAperture, aperture_photometry
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

'''Authors:
		Patrick Tamburo, Boston University, June 2020
	Purpose:
        Chooses suitable reference stars for a target.
	Inputs:
        target (str): The target's full 2MASS name
        source_detect_image_ind (int, optional): The index of the reduced file to use for source detection. By default, will use the first reduced file in the directory.
        radius_check (float, optional): The radius (in pixels) of an aperture that will be placed around potential reference stars to see if there are non-linear pixels.
        non_linear_limit (float, optional): If a potential reference star has a pixel within an aperture of radius radius_check with a value greater than non_linear_limit, it will not be used.
            NOTE: Mimir data becomes non-linear around 4000 ADU.
        dimness_tolerance (float, optional): The minimum fraction of the target's brightness that a reference star may have. 
        closeness_tolerance (float, optional): The closest a potential reference star can be to another source before it's discarded (in pixels). 
        distance_from_target (float, optional): The furthest a potential reference star can be from the target before it's discarded (in pixels).
        lower_left (bool, optional): Whether or not to automatically discard potential references in the lower left quadrant, which is necessary if Mimir's bars problem exists in the data. 
        restore (bool, optional): Whether or not to restore ref_star_chooser output that already exists. 

    Outputs:
		ref stars...
	TODO:
'''

def ref_star_chooser(target, source_detect_image_ind=0, radius_check=6., non_linear_limit=3300., dimness_tolerance=1.0, closeness_tolerance=12., distance_from_target=700., lower_left=False, restore=False):
    
    #Get the target's 'short name'
    short_name = short_name_creator(target)

    #Get your local PINES directory
    pines_path = pines_dir_check()

    if restore:
        output_df = pd.read_csv(pines_path/('Objects/'+short_name+'/sources/target_and_references.csv'))
        print('')
        print('Restoring ref_star_chooser output from {}.'.format(pines_path/('Objects/'+short_name+'/sources/target_and_references.csv')))
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
    phot_dir = pines_path/('Objects/'+short_name+'/sources/')
    sources = detect_sources(reduced_files[source_detect_image_ind], title=reduced_files[source_detect_image_ind].name, plot=False)
    
    #Save detected sources to csv file. 
    print('')
    print('Saving detected sources to {}'.format(phot_dir/(reduced_files[source_detect_image_ind].name.split('.fits')[0]+'_source_detection.csv')))
    #sources.to_csv(phot_dir/(reduced_files[source_detect_image_ind].name.split('.fits')[0]+'_source_detection.csv'))
    with open(phot_dir/(reduced_files[source_detect_image_ind].name.split('.fits')[0]+'_source_detection.csv'),'w') as ofile:
        np.savetxt(ofile, sources.values, fmt='%-4d,%25s,%25s,%25s', header='id,         xcentroid,                      ycentroid,            aperture_sum')
    #print('Saving detected source image to {}'.format(phot_dir/(reduced_files[source_detect_image_ind].name.split('.fits')[0]+'_source_detection.png')))
    print('')
    #plt.savefig(phot_dir/(reduced_files[source_detect_image_ind].name.split('.fits')[0]+'_source_detection.png'))

    target_id = target_finder(sources)
    #ax.plot(sources['xcenter'][target_id], sources['ycenter'][target_id], 'mx')

    #Now determine ids of suitable reference stars.
    image = fits.open(reduced_files[source_detect_image_ind])[0].data
    suitable_ref_ids = []
    bad_ref_ids = [target_id]
    #Get a value for the brightness of the target with a radius_check aperture. We don't want reference stars dimmer than dimness_tolerance * targ_flux_estimates. 
    target_ap = CircularAperture((sources['xcenter'][target_id], sources['ycenter'][target_id]), r=radius_check)
    targ_flux_estimate = aperture_photometry(image, target_ap)['aperture_sum']
    for i in range(len(sources)):
        if i != target_id:
            potential_ref_loc = (sources['xcenter'][i], sources['ycenter'][i])
            #Identify brightnesses of pixels surrounding the potential reference within a 10 pixel radius. 
            ap = CircularAperture(potential_ref_loc, r=radius_check)
            pixels = ap.to_mask()
            sub_image = pixels.cutout(image)

            dists = []
            for j in range(len(sources)):
                dists.append(np.sqrt((sources['xcenter'][i]-sources['xcenter'][j])**2+(sources['ycenter'][i]-sources['ycenter'][j])**2))
            dists = np.array(dists)
            
            non_linear_flag = len(np.where(sub_image > non_linear_limit)[0]) == 0 #Check if potential ref is near non-linear regime.
            dimness_flag = (aperture_photometry(image, ap)['aperture_sum'] > dimness_tolerance*targ_flux_estimate)[0] #Check if potential ref is bright enough. 
            closeness_flag = (len(np.where(dists[np.where(dists != 0)] < closeness_tolerance)[0]) == 0)
            proximity_flag = (np.sqrt((sources['xcenter'][target_id]-sources['xcenter'][i])**2+(sources['ycenter'][target_id]-sources['ycenter'][i])**2) < distance_from_target)
            if lower_left:
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
    #suitable_refs.reset_index()

    output_df = pd.DataFrame(columns = ['Name','Source Detect X','Source Detect Y'])
    output_df = output_df.append({'Name': short_name, 'Source Detect X': sources['xcenter'][target_id], 'Source Detect Y': sources['ycenter'][target_id]}, ignore_index=True)
    for i in range(len(suitable_refs)):
        output_df = output_df.append({'Name': 'Reference '+str(i+1), 'Source Detect X': sources['xcenter'][suitable_ref_ids[i]], 'Source Detect Y': sources['ycenter'][suitable_ref_ids[i]]}, ignore_index=True)
    
    print('Saving target/reference info to {}.'.format(phot_dir/'target_and_references.csv'))
    output_df.to_csv(phot_dir/'target_and_references.csv',index=False)

    stats = sigma_clipped_stats(image)
    fig, ax = plt.subplots(figsize=(10,9))
    ax.set_aspect('equal')
    im = ax.imshow(image, origin='lower', vmin=stats[1], vmax=stats[1]+7*stats[2], cmap='Greys')
    cax = fig.add_axes([0.87, 0.15, 0.035, 0.7])
    fig.colorbar(im, cax=cax, orientation='vertical', label='ADU')
    ax.set_title(reduced_files[source_detect_image_ind].name)
    ax.plot(output_df['Source Detect X'][0], output_df['Source Detect Y'][0], marker='o', color='m', mew=2, ms=12, ls='', mfc='None', label='Target')
    for i in range(1,len(output_df)):
        ax.plot(output_df['Source Detect X'][i], output_df['Source Detect Y'][i], marker='o', color='b', mew=2, ms=12, ls='', mfc='None', label='Reference')
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.1, 1.1))
    print('')
    print('Saving target and references image to {}.'.format(phot_dir/(reduced_files[source_detect_image_ind].name.split('.fits')[0]+'_target_and_refs.png')))
    plt.savefig(phot_dir/(reduced_files[source_detect_image_ind].name.split('.fits')[0]+'_target_and_refs.png'))
    
    return output_df