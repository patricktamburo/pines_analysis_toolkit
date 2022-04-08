from natsort import natsorted
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pdb
import os
from pathlib import Path
from matplotlib.patches import Polygon

from glob import glob
from photutils import DAOStarFinder, aperture_photometry, CircularAperture, Background2D, MedianBackground

from pines_analysis_toolkit.utils import get_source_names, pines_dir_check, pines_log_reader, short_name_creator, pines_login, quick_plot as qp
from pines_analysis_toolkit.data import get_master_synthetic_image, bg_2d, reduce

from astropy.modeling import models, fitting
from astropy.stats import sigma_clip, sigma_clipped_stats, SigmaClip
from astropy.modeling import models, fitting
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from astropy.io import fits
from astropy.visualization import ImageNormalize, ZScaleInterval, LinearStretch
import astropy.units as u 
from astropy.coordinates import SkyCoord

from scipy import signal
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import scipy.optimize as opt
from scipy.stats import sigmaclip

import pandas as pd
pd.options.mode.chained_assignment = None  # Suppress some useless warnings.


matplotlib.use('TkAgg')
plt.ion()

def detect_sources(image_path, seeing_fwhm, edge_tolerance, thresh=6.0, plot=False):
    """Finds sources in a Mimir image.

    :param image_path: path to the image
    :type image_path: pathlib.PosixPath
    :param seeing_fwhm: seeing FWHM in arcsec
    :type seeing_fwhm: float
    :param edge_tolerance: how close a source can be to the edge (pixels)
    :type edge_tolerance: float
    :param thresh: significance threshold, defaults to 6.0
    :type thresh: float, optional
    :param plot: whether or not to plot detected sources, defaults to False
    :type plot: bool, optional
    :return: dataframe of sources
    :rtype: pandas DataFrame
    """

    fwhm = seeing_fwhm/0.579  # FIXED

    # Radius of aperture in pixels for doing quick photometry on detected sources.
    ap_rad = 7

    # Read in the image.
    image = fits.open(image_path)[0].data
    header = fits.open(image_path)[0].header

    # Interpolate nans if any in image.
    kernel = Gaussian2DKernel(x_stddev=0.5)
    image = interpolate_replace_nans(image, kernel)

    # Do a quick background model subtraction (makes source detection easier).
    image = bg_2d(image, box_size=32)

    # Get the sigma_clipped_stats for the image.
    avg, med, std = sigma_clipped_stats(image)

    norm = ImageNormalize(
        data=image, interval=ZScaleInterval(), stretch=LinearStretch())

    if plot:
        title = image_path.name
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 9))
        ax.set_aspect('equal')
        im = ax.imshow(image, origin='lower', norm=norm)
        cax = fig.add_axes([0.94, 0.15, 0.015, 0.7])
        fig.colorbar(im, cax=cax, orientation='vertical', label='ADU')
        ax.set_title(title+' Initial Source Detection')

    print('')
    print('Finding sources in {}.'.format(image_path.name))

    # Detect sources using DAOStarFinder.
    daofind = DAOStarFinder(fwhm=fwhm, threshold=thresh*std, sharplo=0.2)

    initial_sources = daofind(image - med)

    # Do a cut based on source sharpness to get rid of some false detections.
    initial_sources.sort('sharpness')
    bad_sharpness_locs = np.where(initial_sources['sharpness'] < 0.3)[0]
    initial_sources.remove_rows(bad_sharpness_locs)

    # Cut sources that are found within edge_tolerance pix of the edges.
    bad_x = np.where((initial_sources['xcentroid'] < edge_tolerance) | (
        initial_sources['xcentroid'] > 1023-edge_tolerance))[0]
    initial_sources.remove_rows(bad_x)
    bad_y = np.where((initial_sources['ycentroid'] < edge_tolerance) | (
        initial_sources['ycentroid'] > 1023-edge_tolerance))[0]
    initial_sources.remove_rows(bad_y)

    # Cut sources near y = 512, these are frequently bad.
    bad_512 = np.where((initial_sources['ycentroid'] > 506) & (
        initial_sources['ycentroid'] < 518))
    initial_sources.remove_rows(bad_512)

    # Cut sources near the top and bottom edges to avoid issues with the Mimir "ski jump" feature that can sometimes occur.
    bad_ski_jump = np.where((initial_sources['ycentroid'] < 100) | (
        initial_sources['ycentroid'] > 924))
    initial_sources.remove_rows(bad_ski_jump)

    # Do quick photometry on the remaining sources.
    positions = [(initial_sources['xcentroid'][i], initial_sources['ycentroid'][i])
                 for i in range(len(initial_sources))]
    apertures = CircularAperture(positions, r=ap_rad)
    phot_table = aperture_photometry(image-med, apertures)

    # Cut based on brightness.
    phot_table.sort('aperture_sum')
    cutoff = 1*std*np.pi*ap_rad**2
    bad_source_locs = np.where(phot_table['aperture_sum'] < cutoff)
    phot_table.remove_rows(bad_source_locs)
    initial_sources.remove_rows(bad_source_locs)

    if plot:
        # Plot detected sources.
        # TODO: indicate saturated sources, sources near edge, etc. with different color markers.
        ax.plot(phot_table['xcenter'], phot_table['ycenter'],
                'ro', markerfacecolor='none')

    #print('Found {} sources.'.format(len(phot_table)))
    # Resort remaining sources so that the brightest are listed firsts.
    sources = phot_table[::-1].to_pandas()
    return sources

def seeing(log_path):
    """Returns seeing from a PINES observing log. 

    :param log_path: ~/PINES_analysis_toolkit/Logs/ path
    :type log_path: pathlib.PosixPath
    :return: array of seeing values
    :rtype: numpy array 
    """
    df = pines_log_reader(log_path)
    if 'X seeing' in df.keys():
        seeing = np.array(df['X seeing'], dtype=float)
        seeing = np.array(seeing[np.where(~np.isnan(seeing))], dtype=float)
        seeing = seeing[np.where((seeing > 1.2) & (seeing < 7.0))[0]]
        return seeing
        
def seeing_measurer(object_path, box_w=30, plate_scale=0.579, plots=False):
    """Remeasures seeing for all reduced file for the target, and saves a seeing.csv file to the sources directory.

    :param object_path: path to the object directory
    :type object_path: pathlib.PosixPath
    :param box_w: width of boxes to use for measuring source FWHMs, defaults to 30
    :type box_w: int, optional
    :param plate_scale: plate scale of the detector in ''/pixel, defaults to 0.579
    :type plate_scale: float, optional
    :param plots: whether or not to make plots of the fitted gaussians, defaults to False
    :type plots: bool, optional
    """
    def tie_sigma(model):
        return model.x_stddev_1

    red_files = np.array(natsorted(glob(str(object_path/'reduced')+'/*.fits')))

    centroid_df = pines_log_reader(object_path/('sources/target_and_references_centroids.csv'))
    sources = get_source_names(centroid_df)
    centroid_x = np.zeros((len(sources), len(centroid_df)))
    centroid_y = np.zeros((len(sources), len(centroid_df)))
    for i in range(len(sources)):
        centroid_x[i,:] = np.array(centroid_df[sources[i]+' Image X'], dtype='float')
        centroid_y[i,:] = np.array(centroid_df[sources[i]+' Image Y'], dtype='float')
    
    if len(red_files) != len(centroid_x[0,:]):
        print('ERROR: number of reduced files does not equal the number of centroid measurements. Skipping seeing measurements.')
        return

    seeings = np.zeros(len(red_files))

    y, x = np.mgrid[:box_w-2,:box_w-2]

    for i in range(len(red_files)):
        print('{} of {}'.format(i+1, len(red_files)))
        im = fits.open(red_files[i])[0].data
        fwhms = []
        for j in range(len(sources)):
            if np.isnan(centroid_x[j,i]):
                seeings[i] = np.nan
                continue
            #Get the cutout from the centroid position
            cutout = im[int(centroid_y[j,i]-box_w/2):int(centroid_y[j,i]+box_w/2), int(centroid_x[j,i]-box_w/2):int(centroid_x[j,i]+box_w/2)]
            #Interpolate any nans in the cutout
            cutout = interpolate_replace_nans(cutout, Gaussian2DKernel(0.25))
            #If any nans persist, interpolate again.
            if sum(sum(np.isnan(cutout))) > 0:
                cutout = interpolate_replace_nans(cutout, Gaussian2DKernel(0.25))

            #Shave 1 pixel of the edges to limit edge effects.
            cutout = cutout[1:cutout.shape[1]-1, 1:cutout.shape[0]-1]
            #Subtract a background estimate.
            cutout = cutout - np.percentile(cutout,5)

            # Fit with constant, bounds, tied x and y sigmas and outlier rejection:
            gaussian_init = models.Const2D(0.0) + models.Gaussian2D(np.max(cutout),cutout.shape[0]/2, cutout.shape[1]/2, 8/2.355,8/2.355,0)
            gaussian_init.x_stddev_1.min = 1.0/2.355
            gaussian_init.x_stddev_1.max = 20.0/2.355
            gaussian_init.y_stddev_1.min = 1.0/2.355
            gaussian_init.y_stddev_1.max = 20.0/2.355
            gaussian_init.y_stddev_1.tied = tie_sigma
            gaussian_init.theta_1.fixed = True
            fit_gauss = fitting.FittingWithOutlierRemoval(fitting.LevMarLSQFitter(),sigma_clip,niter=3,sigma=3.0)
            gain = 8.21 #e per ADU
            read_noise = 19 #ADU
            weights = gain / np.sqrt(np.absolute(cutout)*gain + (read_noise*gain)**2) #1/sigma for each pixel

            try:
                gaussian, mask = fit_gauss(gaussian_init, x, y, cutout, weights)
            except:
                breakpoint()
            fwhm_x = 2.355*gaussian.x_stddev_1.value * plate_scale
            fwhms.append(fwhm_x)

        v,l,h = sigmaclip(fwhms, 2.5, 2.5)

        if np.isnan(np.nanmean(v)) and not np.isnan(centroid_x[j,i]):
            breakpoint()
        seeings[i] = np.nanmean(v)
        
        if plots:
            norm = ImageNormalize(cutout, interval=ZScaleInterval())
            fig, ax = plt.subplots(1,2,figsize=(10,5))
            ax[0].imshow(cutout, origin='lower', norm=norm)
            ax[1].imshow(gaussian(x,y), origin='lower', norm=norm)
            plt.show()


    seeing_df = pd.DataFrame(seeings,columns=['Seeing'])
    output_path = object_path/('sources/seeing.csv')
    seeing_df.to_csv(output_path, index=0)
    return
    
def background(phot_path):
    """Returns target background values from a PINES photometry file

    :param log_path: photometry file path
    :type log_path: pathlib.PosixPath
    :return: array of background values
    :rtype: numpy array 
    """
    df = pines_log_reader(phot_path)
    sources = get_source_names(df)
    if sources[0]+' Background' in df.keys():
        files = np.array(df['Filename'])
        bg = np.array(df[sources[0]+' Background'], dtype=float)
        files = np.array(files[np.where(~np.isnan(bg))])
        bg = np.array(bg[np.where(~np.isnan(bg))], dtype=float)
        
        if len(files) != len(bg):
            return
        
        #Get exptimes 
        exptimes = np.zeros(len(files))
        bands = np.zeros(len(files), dtype='object')
        for i in range(len(files)):
            path = phot_path.parent.parent/('reduced/'+files[i])
            hdr = fits.open(path)[0].header
            exptimes[i] = hdr['EXPTIME']
            bands[i] = hdr['FILTNME2']

        return bg, exptimes, bands

def bad_shift_identifier(target, date, bad_shift_threshold=200.):
    """Identifies bad shift values in PINES observing logs. 

    :param target: name of the target as it is in the log
    :type target: str
    :param date: date of log to check (YYYYMMDD)
    :type date: str
    :param bad_shift_threshold: pixel threshold above which shifts are considered bad, defaults to 200
    :type bad_shift_threshold: float, optional
    """
    pines_path = pines_dir_check()
    log_path = pines_path/('Logs/'+date+'_log.txt')
    log = pines_log_reader(log_path)
    target_inds = np.where(log['Target'] == target)[0]
    x_shifts = np.array(log['X shift'][target_inds])
    y_shifts = np.array(log['Y shift'][target_inds])

    bad_shift_inds = np.where((x_shifts > bad_shift_threshold) | (
        y_shifts > bad_shift_threshold))[0]
    shift_flags = np.zeros(len(target_inds), dtype=int)
    shift_flags[bad_shift_inds] = 1

    return


def log_out_of_order_fixer(log_path, sftp):
    """Fixes logs with out-of-order/duplicate filenames. 

    :param log_path: path to the log
    :type log_path: pathlib.PosixPath
    :param sftp: sftp connection to the PINES server for downloading missing files
    :type sftp: pysftp connection
    """
    def get_download_path():
        """Returns the default downloads path for linux or windows"""
        if os.name == 'nt':
            import winreg
            sub_key = r'SOFTWARE\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders'
            downloads_guid = '{374DE290-123F-4565-9164-39C4925E467B}'
            with winreg.OpenKey(winreg.HKEY_CURRENT_USER, sub_key) as key:
                location = winreg.QueryValueEx(key, downloads_guid)[0]
            return location
        else:
            return os.path.join(os.path.expanduser('~'), 'downloads')
    myfile = open(log_path, 'r')
    original_lines = myfile.readlines()
    myfile.close()

    # Remove any 'test.fits' lines, if they exist.
    lines = []
    for i in range(len(original_lines)):
        if 'test.fits' in original_lines[i].split(',')[0]:
            continue
        else:
            lines.append(original_lines[i])

    # Sort everything based on file number.
    log_filenums = []
    for i in range(len(lines)):
        if lines[i].strip()[0] == '#':
            log_filenums.append(0)
        else:
            log_filenums.append(int(lines[i].split(',')[0].split('.')[1]))

    sort_inds = np.argsort(log_filenums)
    lines = list(np.array(lines)[sort_inds])
    log_filenums = np.array(log_filenums)[sort_inds]

    num_del = 0
    # Check for lack of entries, multiple entries for this file number in the log.

    for i in range(0, np.max(log_filenums)):
        file_num = i + 1
        entry_inds = np.where(log_filenums == file_num)[0]
        num_entries = len(entry_inds)

        # Check if there are no entries for this file number...if not, create one by finding the file on the PINES server,
        # downloading it, and reading its header.
        if num_entries == 0:
            found = False

            missing_filename = log_path.name.split(
                '_')[0]+'.'+str(file_num).zfill(3)+'.fits'

            if missing_filename != '20200206.607.fits':  # corrupt file
                print('{} missing from log.'.format(missing_filename))
                print('Searching for file on PINES server...')

                pines_raw_path = '/data/raw/mimir/'
                runs = sftp.listdir(pines_raw_path)

                # The file will be in one of these directories if it's on the server
                run_guess = missing_filename[0:6]
                ind = np.where(np.array(runs) == run_guess)[0][0]

                if ind + 1 == len(runs):
                    inds = np.arange(ind-1, ind+1)
                    runs = np.array(runs)[inds]
                    runs = [runs[1], runs[0]]
                else:
                    inds = np.arange(ind-1, ind+2)
                    runs = np.array(runs)[inds]
                    runs = [runs[1], runs[0], runs[2]]

                for jj in range(len(runs)):
                    nights = sftp.listdir(pines_raw_path+'/'+runs[jj])
                    nights = [i for i in nights if i[0] == '2']
                    for kk in range(len(nights)):
                        server_files = sftp.listdir(
                            pines_raw_path+runs[jj]+'/'+nights[kk])
                        if missing_filename in server_files:
                            print('File found: {}'.format(pines_raw_path +
                                  runs[jj]+'/'+nights[kk]+'/'+missing_filename))
                            found = True
                            # Download the file and grab relevant parameters from the header.
                            user_download_path = get_download_path()
                            sftp.get(pines_raw_path+runs[jj]+'/'+nights[kk]+'/' +
                                     missing_filename, user_download_path+'/'+missing_filename)
                            header = fits.open(
                                user_download_path+'/'+missing_filename)[0].header
                            date = header['DATE']
                            if header['OBJECT'] == 'dummy':
                                print(
                                    'ERROR: PINES_watchdog logged "dummy" for object filename; inspect the field and update its name manually.')
                            else:
                                target_name = header['OBJECT'].split(
                                    'J')[0] + ' J' + header['OBJECT'].split('J')[1]
                            filter_name = header['FILTNME2']
                            exptime = str(header['EXPTIME'])
                            airmass = str(header['AIRMASS'])
                            # Use these guesses for shifts/seeings. Will be updated at a later step.
                            x_shift = str(0.0)
                            y_shift = str(0.0)
                            x_seeing = str(2.5)
                            y_seeing = str(2.5)

                            # Add the line to the log.
                            line = pines_logging(missing_filename, date, target_name, filter_name,
                                                 exptime, airmass, x_shift, y_shift, x_seeing, y_seeing, '0', '0')
                            lines.insert(i+1, line)

                            # Delete the file from the Downloads folder.
                            os.remove(user_download_path+'/'+missing_filename)
                            break

            if not found:
                print('{} not found in any raw directories on the PINES server... inspect manually.'.format(
                    missing_filename))
                pdb.set_trace()

        # Check for multiple entries in the log...if so, delete the extras.
        elif num_entries > 1:
            # Only retain the FIRST entry in the log. Assume later entries are mistakes.
            num_del += 1
            for j in range(len(entry_inds)-1, 0, -1):
                del lines[entry_inds[j]]

        # Update the log on disk.
        with open(log_path, 'w') as f:
            for line in lines:
                f.write(line)

    return


def log_updater(date, sftp, shift_tolerance=30., upload=False, force_output_path=''):
    """Updates x_shift and y_shift measurements from a PINES log. These shifts are measured using *full* resolution images, while at the telescope,
        we use *half* resolution images (to save time between exposures). By measuring on full-res images, we get more accurate shifts, which allows 
        us to determine centroids more easily.

    :param date: date of the log whose shifts you want to update (YYYYMMDD)
    :type date: str
    :param sftp: sftp connection to the PINES server
    :type sftp: pysftp connection
    :param shift_tolerance: maximum pixel distance an x/y shift can be before shifts flagged as poor quality, defaults to 30
    :type shift_tolerance: float, optional
    :param upload: [description], whether or not to push the updated log to the PINES server
    :type upload: bool, optional
    :param force_output_path: user-chosen path if you do not want to use the default ~/Documents/PINES_analysis_toolkit/ directory for analysis, defaults to ''
    :type force_output_path: str, optional
    """

    def tie_sigma(model):
        return model.x_stddev_1

    def guide_star_seeing(subframe):
        # subframe = subframe - np.median(subframe)
        subframe = subframe - np.percentile(subframe, 5)
        sub_frame_l = int(np.shape(subframe)[0])
        y, x = np.mgrid[:sub_frame_l, :sub_frame_l]

        # Fit with constant, bounds, tied x and y sigmas and outlier rejection:
        gaussian_init = models.Const2D(0.0) + models.Gaussian2D(subframe[int(sub_frame_l/2), int(
            sub_frame_l/2)], int(sub_frame_l/2), int(sub_frame_l/2), 8/2.355, 8/2.355, 0)
        gaussian_init.x_stddev_1.min = 1.0/2.355
        gaussian_init.x_stddev_1.max = 20.0/2.355
        gaussian_init.y_stddev_1.min = 1.0/2.355
        gaussian_init.y_stddev_1.max = 20.0/2.355
        gaussian_init.y_stddev_1.tied = tie_sigma
        gaussian_init.theta_1.fixed = True
        fit_gauss = fitting.FittingWithOutlierRemoval(
            fitting.LevMarLSQFitter(), sigma_clip, niter=3, sigma=3.0)
        # gaussian, mask = fit_gauss(gaussian_init, x, y, subframe)
        gain = 8.21  # e per ADU
        read_noise = 2.43  # ADU
        # 1/sigma for each pixel
        weights = gain / np.sqrt(np.absolute(subframe)
                                 * gain + (read_noise*gain)**2)
        gaussian, mask = fit_gauss(gaussian_init, x, y, subframe, weights)
        fwhm_x = 2.355*gaussian.x_stddev_1.value
        fwhm_y = 2.355*gaussian.y_stddev_1.value

        x_seeing = fwhm_x * 0.579
        y_seeing = fwhm_y * 0.579
        return(x_seeing, y_seeing)

    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()
    
    log_path = pines_path/('Logs/'+date+'_log.txt')

    # Begin by checking filenames, making sure they're in sequential order, and that there is only one entry for each.
    log_out_of_order_fixer(log_path, sftp)

    log = pines_log_reader(log_path)  # Get telescope log shifts.
    myfile = open(log_path, 'r')
    lines = myfile.readlines()
    myfile.close()

    # Now loop over all files in the log, measure shifts in each file and update the line in the log.
    for i in range(len(log)):
        if (log['Target'][i].lower() != 'flat') & (log['Target'][i].lower() != 'skyflat') & (log['Target'][i].lower() != 'supersky') & (log['Target'][i].lower() != 'flat_on') & (log['Target'][i].lower() != 'flat_off') & (log['Target'][i].lower() != 'dark') & (log['Target'][i].lower() != 'bias') & (log['Target'][i].lower() != 'dummy') & (log['Target'][i].lower() != 'linearity') & (log['Post-processing flag'][i] != 1):
            filename = log['Filename'][i].split('.fits')[0]+'_red.fits'
            target = log['Target'][i]
            short_name = short_name_creator(target)
            image_path = pines_path / \
                ('Objects/'+short_name+'/reduced/'+filename)

            # Figure out which file you're looking at and its position in the log.
            log_ind = np.where(log['Filename'] ==
                               filename.split('_')[0]+'.fits')[0][0]

            # Measure the shifts and get positions of targets.
            (measured_x_shift, measured_y_shift, source_x, source_y,
             check_image) = shift_measurer(target, filename, force_output_path=pines_path)

            if (abs(measured_x_shift) > shift_tolerance) or (abs(measured_y_shift) > shift_tolerance):
                print('Shift greater than {} pixels measured for {} in {}.'.format(
                    shift_tolerance, short_name, image_path.name))
                print('Inspect manually.')
                shift_quality_flag = 1
            elif np.isnan(measured_x_shift) or np.isnan(measured_y_shift):
                raise RuntimeError('Found nans for shifts!')
                shift_quality_flag = 1
            else:
                shift_quality_flag = 0

            # Measure the seeing.
            guide_star_cut = np.where((source_x > 50) & (
                source_x < 975) & (source_y > 50) & (source_y < 975))[0]
            if len(guide_star_cut) != 0:
                x_seeing_array = []
                y_seeing_array = []
                for guide_star_ind in guide_star_cut:
                    guide_star_x_int = int(source_x[guide_star_ind])
                    guide_star_y_int = int(source_y[guide_star_ind])
                    guide_star_subframe = check_image[guide_star_y_int -
                                                      15:guide_star_y_int+15, guide_star_x_int-15:guide_star_x_int+15]
                    (x_seeing, y_seeing) = guide_star_seeing(guide_star_subframe)
                    # Cut unrealistic values/saturated stars.
                    if x_seeing > 1.2 and x_seeing < 7.0:
                        x_seeing_array.append(x_seeing)
                        y_seeing_array.append(y_seeing)
                x_seeing = np.nanmedian(x_seeing_array)
                y_seeing = np.nanmedian(y_seeing_array)
            else:
                # Default to the average PINES value if no sources were found for guiding.
                x_seeing = 2.6
                y_seeing = 2.6

            print('Log line {} of {}.'.format(i+1, len(log)))
            print('Measured x shift: {:4.1f}, measured y shift: {:4.1f}'.format(
                measured_x_shift, measured_y_shift))
            print('Measured seeing: {:4.1f}'.format(x_seeing))
            print('')

            # Overwrite the telescope's logged shifts and seeing values with the new measurements.
            log['X shift'][log_ind] = str(np.round(measured_x_shift, 1))
            log['Y shift'][log_ind] = str(np.round(measured_y_shift, 1))
            log['X seeing'][log_ind] = str(np.round(x_seeing, 1))
            log['Y seeing'][log_ind] = str(np.round(y_seeing, 1))

            # Grab entries for log line.
            filename = log['Filename'][log_ind]
            log_date = log['Date'][log_ind]
            target_name = log['Target'][log_ind]
            filter_name = log['Filt.'][log_ind]
            exptime = log['Exptime'][log_ind]
            airmass = log['Airmass'][log_ind]
            x_shift = log['X shift'][log_ind]
            y_shift = log['Y shift'][log_ind]
            x_seeing = log['X seeing'][log_ind]
            y_seeing = log['Y seeing'][log_ind]
            post_processing_flag = 1
            # Generate line of log text following the PINES telescope log format.
            log_text = pines_logging(filename, log_date, target_name, filter_name, exptime, airmass,
                                     x_shift, y_shift, x_seeing, y_seeing, post_processing_flag, shift_quality_flag)

            # Overwrite the line with the new shifts.
            line_ind = log_ind + 1
            lines[line_ind] = log_text

            # Update the log on disk.
            with open(log_path, 'w') as f:
                for line in lines:
                    f.write(line)

        elif (log['Post-processing flag'][i] == 1):
            print('File already post-processed, skipping. {} of {}'.format(i+1, len(log)))
        else:
            print('File not a science target, skipping. {} of {}.'.format(
                i+1, len(log)))

    if upload:
        sftp.chdir('/data/logs/')
        print('Uploading to /data/logs/{}_log.txt.'.format(date))
        sftp.put(log_path, date+'_log.txt')
    return


def master_synthetic_image_creator(target, image_name, seeing=2.5, sigma_above_bg=5.):
    """Creates a master synthetic image for a PINES target by detecting sources in a reduced image of the field.

    :param target: long name of the target
    :type target: str
    :param image_name: the name of the reduced file you want to use to create the master synthetic image
    :type image_name: str
    :param seeing: FWHM seeing value in the image in arcsec, defaults to 2.5
    :type seeing: float, optional
    :param sigma_above_bg: sigma above background used to detect sources, defaults to 5
    :type sigma_above_bg: float, optional
    """

    def mimir_source_finder(image_path, sigma_above_bg, fwhm, exclude_lower_left=False):
        """Find sources in Mimir images."""

        np.seterr(all='ignore')  # Ignore invalids (i.e. divide by zeros)

        # Find stars in the master image.
        avg, med, stddev = sigma_clipped_stats(
            image, sigma=3.0, maxiters=3)  # Previously maxiters = 5!
        daofind = DAOStarFinder(
            fwhm=fwhm, threshold=sigma_above_bg*stddev, sky=med, ratio=0.8)
        new_sources = daofind(image)
        x_centroids = new_sources['xcentroid']
        y_centroids = new_sources['ycentroid']
        sharpness = new_sources['sharpness']
        fluxes = new_sources['flux']
        peaks = new_sources['peak']

        # Cut sources that are found within 20 pix of the edges.
        use_x = np.where((x_centroids > 20) & (x_centroids < 1004))[0]
        x_centroids = x_centroids[use_x]
        y_centroids = y_centroids[use_x]
        sharpness = sharpness[use_x]
        fluxes = fluxes[use_x]
        peaks = peaks[use_x]
        use_y = np.where((y_centroids > 20) & (y_centroids < 1004))[0]
        x_centroids = x_centroids[use_y]
        y_centroids = y_centroids[use_y]
        sharpness = sharpness[use_y]
        fluxes = fluxes[use_y]
        peaks = peaks[use_y]

        # Also cut using sharpness, this seems to eliminate a lot of false detections.
        use_sharp = np.where(sharpness > 0.5)[0]
        x_centroids = x_centroids[use_sharp]
        y_centroids = y_centroids[use_sharp]
        sharpness = sharpness[use_sharp]
        fluxes = fluxes[use_sharp]
        peaks = peaks[use_sharp]

        if exclude_lower_left:
            # Cut sources in the lower left, if bars are present.
            use_ll = np.where((x_centroids > 512) | (y_centroids > 512))
            x_centroids = x_centroids[use_ll]
            y_centroids = y_centroids[use_ll]
            sharpness = sharpness[use_ll]
            fluxes = fluxes[use_ll]
            peaks = peaks[use_ll]

        # Cut targets whose y centroids are near y = 512. These are usually bad.
        use_512 = np.where(np.logical_or(
            (y_centroids < 510), (y_centroids > 514)))[0]
        x_centroids = x_centroids[use_512]
        y_centroids = y_centroids[use_512]
        sharpness = sharpness[use_512]
        fluxes = fluxes[use_512]
        peaks = peaks[use_512]

        # Cut sources with negative/saturated peaks
        use_peaks = np.where((peaks > 30) & (peaks < 3000))[0]
        x_centroids = x_centroids[use_peaks]
        y_centroids = y_centroids[use_peaks]
        sharpness = sharpness[use_peaks]
        fluxes = fluxes[use_peaks]
        peaks = peaks[use_peaks]

        # Do quick photometry on the remaining sources.
        positions = [(x_centroids[i], y_centroids[i])
                     for i in range(len(x_centroids))]
        apertures = CircularAperture(positions, r=4)
        phot_table = aperture_photometry(image-med, apertures)

        # Cut based on brightness.
        phot_table.sort('aperture_sum')
        cutoff = 1*std*np.pi*4**2
        bad_source_locs = np.where(phot_table['aperture_sum'] < cutoff)
        phot_table.remove_rows(bad_source_locs)

        x_centroids = phot_table['xcenter'].value
        y_centroids = phot_table['ycenter'].value

        return(x_centroids, y_centroids)

    def synthetic_image_maker(x_centroids, y_centroids, fwhm):
        # Construct synthetic images from centroid/flux data.
        synthetic_image = np.zeros((1024, 1024))
        sigma = fwhm/2.355
        for i in range(len(x_centroids)):
            # Cut out little boxes around each source and add in Gaussian representations. This saves time.
            int_centroid_x = int(np.round(x_centroids[i]))
            int_centroid_y = int(np.round(y_centroids[i]))
            y_cut, x_cut = np.mgrid[int_centroid_y-10:int_centroid_y +
                                    10, int_centroid_x-10:int_centroid_x+10]
            dist = np.sqrt((x_cut-x_centroids[i])**2+(y_cut-y_centroids[i])**2)
            synthetic_image[y_cut, x_cut] += np.exp(-(
                (dist)**2/(2*sigma**2)+((dist)**2/(2*sigma**2))))
        return(synthetic_image)

    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    master_synthetic_path = pines_path / \
        ('Calibrations/Master Synthetic Images/'+target+'_master_synthetic.fits')
    image_path = pines_path/('Objects/'+short_name+'/reduced/'+image_name)
    plt.ion()

    seeing = float(seeing)
    daostarfinder_fwhm = seeing*2.355/0.579

    # Open the image and calibration files.
    header = fits.open(image_path)[0].header
    image = fits.open(image_path)[0].data

    # Interpolate over bad pixels
    kernel = Gaussian2DKernel(x_stddev=1)
    image = interpolate_replace_nans(image, kernel)

    # Do a simple 2d background model.
    box_size = 32
    sigma_clip = SigmaClip(sigma=3.)
    bkg_estimator = MedianBackground()
    bkg = Background2D(image, (box_size, box_size), filter_size=(
        3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    image = image - bkg.background

    avg, med, std = sigma_clipped_stats(image)

    # Find sources in the image.
    (x_centroids, y_centroids) = mimir_source_finder(
        image, sigma_above_bg=sigma_above_bg, fwhm=daostarfinder_fwhm)

    # Plot the field with detected sources.
    qp(image)
    plt.plot(x_centroids, y_centroids, 'rx')
    for i in range(len(x_centroids)):
        plt.text(x_centroids[i]+8, y_centroids[i] +
                 8, str(i), color='r', fontsize=14)
    plt.title(
        'Inspect to make sure stars were found!\nO for magnification tool, R to reset view')
    plt.tight_layout()
    plt.show()

    print('')
    print('')
    print('')

    # Prompt the user to remove any false detections.
    ids = input('Enter ids of sources to be removed separated by commas (i.e., 4,18,22). If none to remove, hit enter. To break, ctrl + D. ')
    if ids != '':
        ids_to_eliminate = [int(i) for i in ids.split(',')]
        ids = [int(i) for i in np.linspace(
            0, len(x_centroids)-1, len(x_centroids))]
        ids_to_keep = []
        for i in range(len(ids)):
            if ids[i] not in ids_to_eliminate:
                ids_to_keep.append(ids[i])
    else:
        ids_to_keep = [int(i) for i in np.linspace(
            0, len(x_centroids)-1, len(x_centroids))]
    plt.clf()
    plt.imshow(image, origin='lower', vmin=med, vmax=med+5*std)
    plt.plot(x_centroids[ids_to_keep], y_centroids[ids_to_keep], 'rx')
    for i in range(len(x_centroids[ids_to_keep])):
        plt.text(x_centroids[ids_to_keep][i]+8,
                 y_centroids[ids_to_keep][i]+8, str(i), color='r')

    # Create the synthetic image using the accepted sources.
    synthetic_image = synthetic_image_maker(
        x_centroids[ids_to_keep], y_centroids[ids_to_keep], 8)
    plt.figure(figsize=(9, 7))
    plt.imshow(synthetic_image, origin='lower')
    plt.title('Synthetic image')
    plt.show()

    pdb.set_trace()
    print('')
    print('')
    print('')
    # Now write to a master synthetic image.fits file.
    hdu = fits.PrimaryHDU(synthetic_image)
    hdu.writeto(master_synthetic_path, overwrite=True)
    print('Writing master synthetic image to {}/'.format(master_synthetic_path))


def observed_sample_plots(upload=True):
    """Creates observed sample plots for the PINES website

    :param upload: whether or not to upload to the server, defaults to True
    :type upload: bool, optional
    """

    def plot_mwd(RA, Dec, observed_flag, org=0, title='Mollweide projection', projection='mollweide', observed_plot=0):
        ''' 
        Plots targets on the sky in a 'Mollweide' projection.
        RA, Dec are arrays of the same length.
        RA takes values in [0,360), Dec in [-90,90],
        which represent angles in degrees.
        org is the origin of the plot, 0 or a multiple of 30 degrees in [0,360).
        title is the title of the figure.
        projection is the kind of projection: 'mollweide', 'aitoff', 'hammer', 'lambert'
        '''

        x = np.remainder(RA+360-org, 360)  # shift RA values
        ind = x > 180
        x[ind] -= 360    # scale conversion to [-180, 180]
        x = -x    # reverse the scale: East to the left
        x_tick_labels = np.array(
            [150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])  # Label in degrees
        # x_tick_labels = np.array([150,140,130,120,110,100,90,80,70,60,50,40,30,20,10,0,350,340,330,320,310,300,290,280,270,260,250,240,230,220,210]) #FinerLabel in degrees

        x_tick_labels = np.remainder(x_tick_labels+360+org, 360)
        # x_tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])/15 #Label in hours
        # x_tick_labels = np.remainder(x_tick_labels+24+org/15,24)
        x_tick_labels = [int(i) for i in x_tick_labels]
        fig = plt.figure(figsize=(15*.8, 7*.8))
        ax = fig.add_subplot(111, projection=projection)
        # ax.scatter(np.radians(x),np.radians(Dec),color=color,alpha=0.4,zorder=1, label='Targets')  # convert degrees to radians
        for i in range(len(x)):
            if np.array(observed_flag)[i] == 0:
                color = 'k'
            else:
                color = 'k'
                if observed_plot == 1:
                    color = 'g'  # Turn on observed targest plotting.
            ax.scatter(np.radians(x[i]), np.radians(
                Dec[i]), color=color, alpha=0.4, zorder=1, s=25)
        ax.set_yticklabels(
            [str(int(i))+'$^\circ$' for i in np.round(ax.get_yticks()*180/np.pi)], fontsize=15)
        ax.title.set_fontsize(20)
        ax.set_xlabel('RA')
        ax.xaxis.label.set_fontsize(20)
        ax.set_ylabel("Dec")
        ax.yaxis.label.set_fontsize(20)
        # we add the scale on the x axis
        ax.set_xticklabels([], fontsize=16)
        ax.grid(True, alpha=0.3)
        month_texts = ['Sep', 'Aug', 'Jul', 'Jun', 'May',
                       'Apr', 'Mar', 'Feb', 'Jan', 'Dec', 'Nov', 'Oct']
        for i in range(len(month_texts)):
            ax.text(-180*np.pi/180+15*np.pi/180+30*np.pi/180*i, -35*np.pi /
                    180, month_texts[i], ha='center', va='center', fontsize=14)
        for i in range(len(x_tick_labels)):
            ax.text(-150*np.pi/180+30*np.pi/180*i, -22.5*np.pi/180,
                    str(x_tick_labels[i])+'$^\circ$', ha='center', va='center', fontsize=15)

        # Plot monsoon season.
        monsoon_x_vertices = np.array([-150, -150, -90, -90, -150])*np.pi/180
        monsoon_y_vertices = np.array([-90, 90, 90, -90, -90])*np.pi/180
        monsoon_polygon = Polygon(np.array([[monsoon_x_vertices[i], monsoon_y_vertices[i]] for i in range(
            len(monsoon_x_vertices))]), color='r', alpha=0.15, label='Flagstaff monsoon season')
        ax.add_patch(monsoon_polygon)
        plt.show()
        return ax

    '''Plots the current sample as given in 'PINES sample.xlsx' on Google drive and uploads to the PINES website.'''
    pines_path = pines_dir_check()
    sample_path = pines_path/('Misc/PINES Sample.xlsx')
    print('Make sure an up-to-date copy of PINES Sample.xlsx exists in {}.'.format(pines_path/'Misc/'))
    print('Download from the PINES Google Drive.\n')

    df = pd.read_excel(sample_path)
    df = df.dropna(how='all')  # Remove rows that are all NaNs.

    good_locs = np.where(df['Good'] == 1)[0]  # Get only "good" targets
    ra = np.array(df['RA (deg)'][good_locs])
    dec = np.array(df['Dec (deg)'][good_locs])
    group_ids = df['Group ID'][good_locs]
    observed_flag = df['Observed?'][good_locs]
    # Get the groups that have been observed.
    observed_groups = np.unique(
        np.array(group_ids)[np.where(observed_flag != 0)[0]])
    number_observed = len(np.array(group_ids)[np.where(observed_flag != 0)[0]])

    # Plot 1: Sky map of good targets based on group.
    print('Updating sky plot...')
    ax = plot_mwd(ra, dec, observed_flag, org=180,
                  projection='mollweide', observed_plot=1)
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc=1,
              bbox_to_anchor=(1.1, 1.1), fontsize=16)
    ax.grid(alpha=0.2)

    group_id_inds = np.arange(0, max(group_ids)+1)

    # Now loop over group_id inds, and draw boundaries around each group.
    for i in group_id_inds:
        targs_in_group = np.where(group_ids == i)[0]
        try:
            cluster_coords = np.array([[ra[i], dec[i]]
                                      for i in targs_in_group])
        except:
            pdb.set_trace()
        hull = ConvexHull(cluster_coords)
        for s in range(len(hull.simplices)):
            simplex = hull.simplices[s]
            # shift RA values
            x = np.remainder(cluster_coords[simplex, 0]+360-180, 360)
            ind = x > 180
            x[ind] -= 360    # scale conversion to [-180, 180]
            x = -x    # reverse the scale: East to the left
            if i in observed_groups:
                color = 'g'
                ax.plot(x*np.pi/180, cluster_coords[simplex, 1]*np.pi/180,
                        color=color, lw=2, zorder=0, alpha=0.6, label='Observed')
            else:
                color = 'k'
                ax.plot(x*np.pi/180, cluster_coords[simplex, 1]*np.pi/180,
                        color=color, lw=2, zorder=0, alpha=0.6, label='Not yet observed')

    ax.grid(alpha=0.4)
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))

    ax.legend(by_label.values(), by_label.keys(),
              loc=1, bbox_to_anchor=(0.65, 0.225))
    ax.set_title('PINES sample \n '+str(int(max(group_ids)+1)) +
                 ' groups, '+str(len(good_locs))+' targets', fontsize=20)
    plt.tight_layout()
    sky_map_output_path = pines_path/('Misc/updated_sky_plot.png')
    plt.savefig(sky_map_output_path, dpi=300)
    plt.close()

    ntargs = len(df)
    # Now do magnitude/SpT histograms
    print('Updating target histograms...')
    mags = np.zeros(ntargs)
    observed_SpTs = []
    observed_mags = []
    SpT = []
    for i in range(ntargs):
        try:
            #mags[i] = float(df['2MASS H'][i][0:6])
            mags[i] = float(df['2MASS J'][i][0:6])
            SpT.append(df['SpT'][i])
            if df['Observed?'][i] != 0:
                observed_SpTs.append(df['SpT'][i])
                observed_mags.append(mags[i])
        # Some values don't follow the normal +/- convention (they were upper limits in the Gagne sheet), so have to read them in differently.
        except:
            #mags[i] = float(df['2MASS H'][i])
            mags[i] = float(df['2MASS J'][i])
            SpT.append(df['SpT'][i])
            if df['Observed?'][i] != 0:
                observed_SpTs.append(df['SpT'][i])
                observed_mags.append(mags[i])

    mags = mags[good_locs]
    SpT = np.array(SpT)
    observed_SpTs = np.array(observed_SpTs)
    observed_mags = np.array(observed_mags)
    SpT = SpT[good_locs]

    SpT_number = np.zeros(ntargs)
    observed_SpT_numbers = []
    for i in range(ntargs):
        if df['SpT'][i][0] == 'L':
            SpT_number[i] = float(df['SpT'][i][1:])
            if df['Observed?'][i] != 0:
                observed_SpT_numbers.append(SpT_number[i])
        else:
            SpT_number[i] = 10 + float(df['SpT'][i][1:])
            if df['Observed?'][i] != 0:
                observed_SpT_numbers.append(SpT_number[i])
    SpT_number = SpT_number[good_locs]
    SpT_number = np.array(SpT_number)
    observed_SpT_numbers = np.array(observed_SpT_numbers)

    median_mag = np.median(mags)

    scale_factor = 0.5
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(
        18*scale_factor, 15*scale_factor))
    bins = np.array([11.25, 11.75, 12.25, 12.75, 13.25, 13.75,
                    14.25, 14.75, 15.25, 15.75, 16.25, 16.75]) - 0.25
    ax[0].hist(mags, bins=bins, histtype='step',
               lw=3, ls='--', label='Full sample')
    ax[0].hist(observed_mags, bins=bins, histtype='bar',
               label='Observed sample', color='tab:blue')
    ax[0].axvline(median_mag, color='r',
                  label='Median $m_J$ = {:2.1f}'.format(median_mag))
    ticks = [11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5]
    ax[0].plot()
    ax[0].set_xticks(ticks)
    ax[0].set_xticklabels([str(i) for i in ticks])
    ax[0].set_xlabel('$m_J$', fontsize=20)
    ax[0].set_ylabel('Number of targets', fontsize=20)
    ax[0].tick_params(axis='both', which='major', labelsize=16)
    ax[0].legend(fontsize=16, loc='upper left')
    # ax[0].grid(alpha=0.2)

    ax[1].hist(SpT_number, bins=np.arange(-0.5, max(SpT_number)+0.5, 1),
               histtype='step', lw=3, color='orange', ls='--', label='Full sample')
    ax[1].hist(observed_SpT_numbers, bins=np.arange(-0.5, max(SpT_number)+0.5, 1),
               histtype='bar', lw=3, color='orange', label='Observed sample')
    ticks = np.arange(0, max(SpT_number), 1)
    ax[1].set_xticks(ticks)
    ax[1].set_xticklabels(['L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6',
                          'L7', 'L8', 'L9', 'T0', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7'])
    ax[1].set_xlabel('Spectral Type', fontsize=20)
    ax[1].set_ylabel('Number of targets', fontsize=20)
    ax[1].tick_params(axis='both', which='major', labelsize=16)
    ax[1].legend(fontsize=16, loc='upper right')
    # ax[1].grid(alpha=0.2)

    plt.tight_layout()
    histogram_output_path = pines_path/'Misc/target_histograms.png'
    plt.savefig(histogram_output_path, dpi=300)
    plt.close()

    # Edit the observing.html page to update the number of observed targets.
    print('Updating observing.html...')
    if not (pines_path/'Misc/observing.html').exists():
        print('Grabbing copy of observing.html from the PINES server.')
        sftp = pines_login()
        sftp.chdir('/web')
        remote_path = '/web/observing.html'
        local_path = pines_path/('Misc/observing.html')
        sftp.get(remote_path, local_path)
        sftp.close()

    with open(str(pines_path/('Misc/observing.html')), 'r') as f:
        lines = f.readlines()

    edit_line_ind = np.where(
        ['To date, PINES has observed' in i for i in lines])[0][0]
    edit_line = lines[edit_line_ind]
    edit_line = edit_line.replace(edit_line.split(
        '<u>')[1].split('</u>')[0], str(number_observed))
    lines[edit_line_ind] = edit_line
    with open(str(pines_path/('Misc/observing.html')), 'w') as f:
        f.writelines(lines)

    if upload:
        sftp = pines_login()
        print('Uploading plots and observing.html to the PINES server.')
        sftp.chdir('/web/images')
        sftp.put(sky_map_output_path, '/web/images/updated_sky_plot.png')
        sftp.put(histogram_output_path, '/web/images/target_histograms.png')
        sftp.chdir('/web')
        sftp.put(pines_path/('Misc/observing.html'), '/web/observing.html')
        print('PINES website updated!')


def pines_logging(filename, date, target_name, filter_name, exptime, airmass, x_shift, y_shift, x_seeing, y_seeing, post_processing_flag,   shift_quality_flag):
    """Exports a line of log text in the same style as the telescope observing logs. 

    :param filename: string of raw filename (e.g., '20210214.203.fits')
    :type filename: str
    :param date: date of image (YYYYMMDD)
    :type date: str
    :param target_name: taget's long name
    :type target_name: str
    :param filter_name: filter name (e.g., 'J' or 'H')
    :type filter_name: str
    :param exptime: exposure time
    :type exptime: str
    :param airmass: airmass
    :type airmass: str
    :param x_shift: x shift in pixels
    :type x_shift: str
    :param y_shift: y shift in pixels
    :type y_shift: str
    :param x_seeing: x seeing in arcsec
    :type x_seeing: str
    :param y_seeing: y seeing in arcsec
    :type y_seeing: str
    :param post_processing_flag: 0 (not yet post-processed) or 1 (post-processed)
    :type post_processing_flag: int
    :param shift_quality_flag: 0 (good shifts) or 1 (bad shifts)
    :type shift_quality_flag: int
    """
    try:
        log_text = ' {:<19}, {:<20}, {:<30}, {:<6}, {:<8}, {:<8}, {:<8}, {:<8}, {:<9}, {:<7}, {:<21}, {:<20}\n'.format(filename, date, target_name,
                                                                                                                       filter_name, str(
                                                                                                                           exptime),
                                                                                                                       str(airmass), str(
                                                                                                                           x_shift),
                                                                                                                       str(
                                                                                                                           y_shift),
                                                                                                                       str(
                                                                                                                           x_seeing),
                                                                                                                       str(
                                                                                                                           y_seeing),
                                                                                                                       str(
                                                                                                                           post_processing_flag),
                                                                                                                       str(shift_quality_flag))
    except:
        pdb.set_trace()

    return log_text


def shift_measurer(target, image_name, num_sources=15, closeness_tolerance=10., force_output_path=''):
    """Measure shifts between an image and the master synthetic image for a target. 

    :param target: target's full 2MASS name 
    :type target: str
    :param image_name: name of the image whose shifts you want to measure
    :type image_name: str
    :param num_sources: number of sources to use for measuring shift, defaults to 15
    :type num_sources: int, optional
    :param closeness_tolerance: closest targets can be in pixels and still be considered for correlation sources, defaults to 10.
    :type closeness_tolerance: float, optional
    :param force_output_path: user-chosen path if you do not want to use the default ~/Documents/PINES_analysis_toolkit/ directory for analysis, defaults to ''
    :type force_output_path: str, optional
    :return: x/y shift, source x/y, check_image
    :rtype: [type]
    """

    def corr_shift_determination(corr):
        # Measure shift between the check and master images by fitting a 2D gaussian to corr. This gives sub-pixel accuracy.
        # Find the pixel with highest correlation, then use this as estimate for gaussian fit.
        y_max, x_max = np.unravel_index(np.argmax(corr), corr.shape)
        y, x = np.mgrid[y_max-10:y_max+10, x_max-10:x_max+10]
        corr_cut = corr[y, x]
        gaussian_init = models.Gaussian2D(
            np.max(corr_cut), x_max, y_max, 8/2.355, 8/2.355, 0)
        fit_gauss = fitting.LevMarLSQFitter()
        gaussian = fit_gauss(gaussian_init, x, y, corr_cut)
        fit_x = gaussian.x_mean.value
        fit_y = gaussian.y_mean.value

        x_shift = fit_x - 1024
        y_shift = fit_y - 1024
        return(x_shift, y_shift)

    #image_name = '20201005.354_red.fits'

    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()    
    
    short_name = short_name_creator(target)
    synthetic_filename = target.replace(' ', '')+'_master_synthetic.fits'
    synthetic_path = pines_path / \
        ('Calibrations/Master Synthetic Images/'+synthetic_filename)
    
    image_path = pines_path/('Objects/'+short_name+'/reduced/'+image_name)
    # Check for the appropriate master synthetic image on disk. If it's not there, get from PINES server.
    if not synthetic_path.exists():
        print('No master image in {} for {}.'.format(pines_path, short_name))
        print('Download from PINES server.')
        sftp = pines_login()
        get_master_synthetic_image(sftp, target)

    master_synthetic_image = fits.open(synthetic_path)[0].data
    
    # If the reduced image doesn't exist, download raw from PINES server and reduce it.
    if not image_path.exists():
        missing_filename = image_path.name.split('_')[0]+'.fits'
        try:
            sftp
        except NameError:
            print('Connect to PINES server to download missing image')
            sftp = pines_login()
        sftp.chdir('/data/raw/mimir/')
        runs = sftp.listdir()
        run_guess = missing_filename[0:6]
        ind = np.where(np.array(runs) == run_guess)[0][0]

        if ind + 1 == len(runs):
            inds = np.arange(ind-1, ind+1)
            runs = np.array(runs)[inds]
            runs = [runs[1], runs[0]]
        else:
            inds = np.arange(ind-1, ind+2)
            runs = np.array(runs)[inds]
            runs = [runs[1], runs[0], runs[2]]

        for jj in range(len(runs)):
            nights = sftp.listdir('/data/raw/mimir/'+runs[jj])
            nights = [i for i in nights if i[0] == '2']
            for kk in range(len(nights)):
                server_files = sftp.listdir(
                    '/data/raw/mimir/'+runs[jj]+'/'+nights[kk])
                if missing_filename in server_files:
                    print('File found: {}'.format('/data/raw/mimir/' +
                          runs[jj]+'/'+nights[kk]+'/'+missing_filename))
                    found = True
                    # Download the file and grab relevant parameters from the header.
                    sftp.get('/data/raw/mimir/'+runs[jj]+'/'+nights[kk]+'/'+missing_filename, pines_path/('Objects/'+short_name+'/raw/'+missing_filename))
                    break
        reduce(short_name)

    # Read in the check image and interpolate/background subtract to make source detection easier.
    check_image = fits.open(image_path)[0].data
    check_image = interpolate_replace_nans(
        check_image, kernel=Gaussian2DKernel(x_stddev=0.5))
    bg_2d = Background2D(check_image, box_size=64)
    check_image = check_image - bg_2d.background

    # Read in the log and figure out what seeing FWHM to use.
    date = image_name.split('.')[0]
    log = pines_log_reader(pines_path/('Logs/'+date+'_log.txt'))
    raw_filename = image_name.split('_')[0]+'.fits'
    ind = np.where(log['Filename'] == raw_filename)[0][0]
    seeing = float(log['X seeing'][ind])

    if (seeing <= 1.0) or (seeing >= 7.0) or (np.isnan(seeing)):
        if ind >= 5:
            seeing = np.nanmedian(
                np.array(log['X seeing'][ind-5:ind], dtype='float'))
        else:
            seeing = 2.6

    # Find sources in the image.
    sources = detect_sources(image_path, seeing, edge_tolerance=10, thresh=3.5)

    # Comb the returned sources and cut any that are too close to another one.
    # This can happen if the actual seeing differs from that recorded in the log.
    bad_inds = []
    for i in range(len(sources)):
        if i not in bad_inds:
            x = sources['xcenter'][i]
            y = sources['ycenter'][i]
            dists = np.array(
                np.sqrt((sources['xcenter']-x)**2 + (sources['ycenter']-y)**2))
            duplicate_inds = np.where(
                (dists < closeness_tolerance) & (dists != 0))[0]
            if len(duplicate_inds) >= 1:
                bad_inds.extend(duplicate_inds)

    ap_sum = np.array(sources['aperture_sum'])
    source_x = np.array(sources['xcenter'])
    source_y = np.array(sources['ycenter'])

    ap_sum = np.delete(ap_sum, bad_inds)
    source_x = np.delete(source_x, bad_inds)
    source_y = np.delete(source_y, bad_inds)

    sort_inds = np.argsort(ap_sum)[::-1]
    source_x = source_x[sort_inds[0:num_sources]]
    source_y = source_y[sort_inds[0:num_sources]]

    # if raw_filename == '20201005.354.fits':
    #     qp(check_image)
    #     plt.plot(source_x, source_y, 'rx')
    #     pdb.set_trace()

    check_synthetic_image = synthetic_image_maker(source_x, source_y)

    # Measure the shift between the synthetic images.
    corr = signal.fftconvolve(master_synthetic_image,
                              check_synthetic_image[::-1, ::-1])

    (x_shift, y_shift) = corr_shift_determination(corr)
    #print('(X shift, Y shift): ({:3.1f}, {:3.1f})'.format(x_shift, -y_shift))
    # print('')

    return (x_shift, -y_shift, source_x, source_y, check_image)


def synthetic_image_maker(x_centroids, y_centroids, fwhm=2.5):
    """Construct a synthetic image from centroid data

    :param x_centroids: array of x positions
    :type x_centroids: numpy array 
    :param y_centroids: array of y positions
    :type y_centroids: numpy array
    :param fwhm: fwhm of synthetic sources, defaults to 8
    :type fwhm: float, optional
    """
    # Construct synthetic images from centroid/flux data.
    synthetic_image = np.zeros((1024, 1024))
    sigma = fwhm/2.355
    for i in range(len(x_centroids)):
        # Cut out little boxes around each source and add in Gaussian representations. This saves time.
        int_centroid_x = int(np.round(x_centroids[i]))
        int_centroid_y = int(np.round(y_centroids[i]))
        y_cut, x_cut = np.mgrid[int_centroid_y-10:int_centroid_y +
                                10, int_centroid_x-10:int_centroid_x+10]
        dist = np.sqrt((x_cut-x_centroids[i])**2+(y_cut-y_centroids[i])**2)
        synthetic_image[y_cut,
                        x_cut] += np.exp(-((dist)**2/(2*sigma**2)+((dist)**2/(2*sigma**2))))
    return(synthetic_image)

def group_separation_measurer(group_id):
    """Calculates the separation between group members in degrees.

    :param group_id: The group ID from the PINES sample
    :type group_id: float
    """

    pines_path = pines_dir_check()
    sample_path = pines_path/('Misc/PINES Sample.csv')
    if not os.path.exists(sample_path):
        print('Download a csv of the PINES sample from Google Drive.')
        return 
    
    sample_df = pd.read_csv(sample_path, converters={'RA (h:m:s)':str, 'Dec (d:m:s)':str})
    group_inds = np.where((sample_df['Group ID'] == group_id) & (sample_df['Good'] == 1))[0]
    n_targs = len(group_inds)
    ra_array = np.array(sample_df['RA (h:m:s)'][group_inds])
    dec_array = np.array(sample_df['Dec (d:m:s)'][group_inds])

    for i in range(len(group_inds)):
        s1 = SkyCoord(ra_array[i]+' '+dec_array[i], unit=(u.hourangle, u.deg))
        s2 = SkyCoord(ra_array[(i+1)%n_targs]+' '+dec_array[(i+1)%n_targs], unit=(u.hourangle, u.deg))
        sep = s1.separation(s2)
        print('Object {} separated from object {} by {:1.1f}.'.format(i+1, (i+2)%(n_targs+1), sep))
    return