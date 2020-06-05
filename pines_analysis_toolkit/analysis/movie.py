from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import pdb
import glob 
from astropy.stats import sigma_clipped_stats
import pickle
from photutils import CircularAperture
from photutils import DAOStarFinder
import warnings

def movie(target_name,run,ap_rad,dates_to_examine=[],target_number=0,box_size=300,show_bpm=0,show_centroid=1):
    '''
        Author: Patrick Tamburo, BU, Jan 2020
        Purpose: Shows a movie of a cutout around a target/reference of your choosing. Useful for seeing if the source falls on bad pixels over the time series, and if your centroids are good.  
        Inputs: 
            target_name (str): The short name of your target (i.e. '2MASS 0046+0715')
            local_path (str): The top-level path to your directory for this target. 
            ap_rad (float): The aperture radius to overplot (useful for seeing if you used too small/too large of an aperture).
            dates_to_examine (list of strings): The dates you want to view a movie of. By default, empty (and will do all dates!). I.e., ['20191118','20191119']
            target_number (int): The id of the source you want a movie of. 0 = target (default), all others correspond to reference stars. 
            box_size (int): The size (in pixels) of the cutout you want.  
            show_bpm (bool): Whether or not to overplot the bad pixel mask.
            show_centroid(bool): Whether or not to overplot the measured centroid positions.
        #TODO: 
            Make it so it can save as a gif? 
    '''
    local_path = '/Users/tamburo/Documents/Data/PINES/Objects/'+target_name+'/'
    warnings.filterwarnings('ignore', category=UserWarning, append=True) #This turns of NaN warnings in sigma_clipped_stats, otherwise we'd get a warning every line. 
    reduced_path = local_path+'domeflat_reduced/'
    reduced_files = sorted(glob.glob(reduced_path+'*red.fits'))
    aper_phot_path = local_path+'/aper_phot/'
    centroid_path = aper_phot_path+'centroids/'
    positions = pickle.load(open(centroid_path+target_name+'_centroids.p','rb'))
    x_positions = positions['X pos'][target_number]
    y_positions = positions['Y pos'][target_number]
    initial_position = (int(x_positions[0]),int(y_positions[0]))

    #Get list of files
    all_frame_list = np.array([])
    for i in range(np.size(reduced_files)):
        file_name = reduced_files[i].split('/')[-1]
        all_frame_list = np.append(all_frame_list, file_name)
    num_files = len(reduced_files)

    dates = np.array([reduced_files[i].split('/')[-1].split('.')[0] for i in range(len(reduced_files))])

    bpm = (1-pickle.load(open('/Users/tamburo/Documents/Data/PINES/Calibrations/Bad Pixel Masks/bpm.p','rb'))).astype('bool') 
    bad_x = np.where(bpm == 0)[1]
    bad_y = np.where(bpm == 0)[0]

    plt.ion()
    fig, ax = plt.subplots(1,1,figsize=(8,7))
    for i in range(len(reduced_files)):

        if len(dates_to_examine) != 0:
            if reduced_files[i].split('/')[-1].split('.')[0] in dates_to_examine:
                image_path = reduced_files[i]
                title = image_path.split('/')[-1]
                image = fits.open(image_path)[0].data
                header = fits.open(image_path)[0].header
                frame = image[initial_position[1]-int(box_size/2):initial_position[1]+int(box_size/2),initial_position[0]-int(box_size/2):initial_position[0]+int(box_size/2)]
                ap = CircularAperture((x_positions[i],y_positions[i]),r=ap_rad)
                avg,med,std = sigma_clipped_stats(image)
                im = ax.imshow(image,origin='lower',vmin=med,vmax=med+5*std)
                if show_bpm:
                    ax.plot(bad_x,bad_y,'rx',alpha=0.7)
                ax.set_xlim(initial_position[0]-int(box_size/2),initial_position[0]+int(box_size/2))
                ax.set_ylim(initial_position[1]-int(box_size/2),initial_position[1]+int(box_size/2))
                cb = fig.colorbar(im, orientation='vertical',label='Counts')
                if show_centroid:
                    ap.plot(color='b',lw=2)
                    ax.plot(x_positions[i],y_positions[i],'bx')
                ax.set_title(title)
                plt.pause(0.01)
                ax.cla()
                cb.remove()
        else:
            image_path = reduced_files[i]
            title = image_path.split('/')[-1]
            image = fits.open(image_path)[0].data
            header = fits.open(image_path)[0].header
            frame = image[initial_position[1]-int(box_size/2):initial_position[1]+int(box_size/2),initial_position[0]-int(box_size/2):initial_position[0]+int(box_size/2)]
            ap = CircularAperture((x_positions[i],y_positions[i]),r=ap_rad)
            avg,med,std = sigma_clipped_stats(image)
            im = ax.imshow(image,origin='lower',vmin=med,vmax=med+5*std)
            if show_bpm:
                ax.plot(bad_x,bad_y,'rx',alpha=0.7)
            ax.set_xlim(initial_position[0]-int(box_size/2),initial_position[0]+int(box_size/2))
            ax.set_ylim(initial_position[1]-int(box_size/2),initial_position[1]+int(box_size/2))
            cb = fig.colorbar(im, orientation='vertical',label='Counts')
            if show_centroid:
                ap.plot(color='b',lw=2)
                ax.plot(x_positions[i],y_positions[i],'bx')
            ax.set_title(title)
            plt.pause(0.1)
            ax.cla()
            cb.remove()
