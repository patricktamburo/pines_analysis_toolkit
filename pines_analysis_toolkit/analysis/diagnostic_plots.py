import matplotlib.pyplot as plt 
import pdb 
import numpy as np 
from natsort import natsorted
from pines_analysis_toolkit.utils.get_source_names import get_source_names
from pines_analysis_toolkit.analysis.night_splitter import night_splitter
from pines_analysis_toolkit.analysis.block_splitter import block_splitter
from pines_analysis_toolkit.analysis.analysis_plots import standard_x_range
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from scipy.interpolate import CubicSpline
import os 
from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader

def plot_style():
    title_size = 20 
    axis_title_size = 20
    axis_ticks_font_size = 16
    legend_font_size = 12
    return (title_size, axis_title_size, axis_ticks_font_size, legend_font_size)

def seeing_plot(target, centroided_sources):
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    #Get plot style parameters. 
    title_size, axis_title_size, axis_ticks_font_size, legend_font_size = plot_style()
    
    #Get list of souce names in the centroid output.
    source_names = get_source_names(centroided_sources)
    centroided_sources.columns = centroided_sources.keys().str.strip()

    #Read in times and seeing values. 
    times_full = np.array(centroided_sources['Time (JD UTC)'])
    seeing = np.array(centroided_sources['Seeing'])

    #Split up times by nights.
    night_inds = night_splitter(times_full)
    num_nights = len(night_inds)
    times_nights = [times_full[night_inds[i]] for i in range(num_nights)]
    standard_x = standard_x_range(times_nights)

    fig, ax = plt.subplots(nrows=1, ncols=num_nights, figsize=(17,5), sharey=True)
    for i in range(num_nights):
        if i == 0:
            ax[i].set_ylabel('Seeing (")', fontsize=axis_title_size)

        inds = night_inds[i]
        ax[i].plot(times_nights[i], seeing[inds], marker='.', linestyle='', alpha=0.3, label='Unbinned seeing')
        ax[i].tick_params(labelsize=axis_ticks_font_size)
        ax[i].set_xlabel('Time (JD UTC)', fontsize=axis_title_size)
        ax[i].grid(alpha=0.2)
        ax[i].set_xlim(np.mean(times_nights[i])-standard_x/2, np.mean(times_nights[i])+standard_x/2)

        #bin
        block_inds = block_splitter(times_nights[i])
        block_x = np.zeros(len(block_inds))
        block_y = np.zeros(len(block_inds))
        block_y_err = np.zeros(len(block_inds))
        for j in range(len(block_inds)):
            block_x[j] = np.mean(times_nights[i][block_inds[j]])
            block_y[j] = np.mean(seeing[inds][block_inds[j]])
            block_y_err[j] = np.std(seeing[inds][block_inds[j]]) / np.sqrt(len(seeing[inds][block_inds[j]]))

        ax[i].errorbar(block_x, block_y, block_y_err, marker='o', linestyle='', color='tab:blue', ms=8, mfc='none', mew=2, label='Binned seeing')
        
        #Interpolate each night's seeing.
        fit_times = np.linspace(block_x[0], block_x[-1], 1000)
        interp = CubicSpline(block_x, block_y)
        interp_fit = interp(fit_times)
        ax[i].plot(fit_times, interp_fit, color='r', lw=2, zorder=0, alpha=0.7, label='CS Interpolation')

    ax[i].legend(bbox_to_anchor=(1.01, 0.5), fontsize=legend_font_size)
    plt.suptitle(short_name+' Seeing Measurements', fontsize=title_size)
    plt.subplots_adjust(left=0.07, hspace=0.05, wspace=0.05, top=0.92, bottom=0.17, right=0.85)

    output_filename = pines_path/('Objects/'+short_name+'/analysis/diagnostic_plots/'+short_name+'_seeing.png')
    plt.savefig(output_filename, dpi=300)
    return

def relative_cutout_position_plot(target, centroided_sources):
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)

    #Get plot style parameters. 
    title_size, axis_title_size, axis_ticks_font_size, legend_font_size = plot_style()
    
    #Get list of souce names in the centroid output.
    source_names = get_source_names(centroided_sources)
    centroided_sources.columns = centroided_sources.keys().str.strip()

    #Get times from the centroid output and split them by night. 
    times_full = np.array(centroided_sources['Time (JD UTC)'])
    night_inds = night_splitter(times_full)
    num_nights = len(night_inds)
    times_nights = [times_full[night_inds[i]] for i in range(num_nights)]
    standard_x = standard_x_range(times_nights)

    #Get the box size (I don't like that this is being determined by using the mean of the data...output it from centroider?)    
    box_w = int(np.round(2*np.mean(np.array(centroided_sources['Reference 1 Cutout X'])),0))

    #X plot
    fig, ax = plt.subplots(nrows=2, ncols=num_nights, figsize=(17,9), sharey=True)
    markers = ['+', 'x', '*', 'X']
    for j in range(num_nights):
        inds = night_inds[j]
        if j == 0: 
            ax[0,j].set_ylabel('Cutout X Position', fontsize=axis_title_size)
            ax[1,j].set_ylabel('Cutout Y Position', fontsize=axis_title_size)

        for i in range(len(source_names)):
            cutout_x = np.array(centroided_sources[source_names[i]+' Cutout X'][inds])
            cutout_y = np.array(centroided_sources[source_names[i]+' Cutout Y'][inds])

            if i == 0:
                marker = 'o'
            else:
                marker = markers[(i-1) % len(markers)]

            ax[0,j].plot(times_nights[j], cutout_x, marker=marker, label=source_names[i], linestyle='')
            ax[1,j].plot(times_nights[j], cutout_y, marker=marker, linestyle='')

        ax[0,j].tick_params(labelsize=axis_ticks_font_size)
        ax[0,j].set_xticklabels([])
        ax[0,j].axhline(box_w/2, zorder=0, color='r', label='Central cutout pixel', lw=2)
        ax[0,j].set_xlim(np.mean(times_nights[j])-standard_x/2, np.mean(times_nights[j])+standard_x/2)
        ax[0,j].grid(alpha=0.2)
        ax[1,j].tick_params(labelsize=axis_ticks_font_size)
        ax[1,j].axhline(box_w/2, zorder=0, color='r', label='Central cutout pixel', lw=2)
        ax[1,j].set_xlim(np.mean(times_nights[j])-standard_x/2, np.mean(times_nights[j])+standard_x/2)
        ax[1,j].set_xlabel('Time (JD UTC)', fontsize=axis_title_size)
        ax[1,j].grid(alpha=0.2)

    plt.suptitle(short_name+' Cutout Centroid Positions', fontsize=title_size)
    plt.subplots_adjust(left=0.07, hspace=0.05, wspace=0.05, top=0.92, bottom=0.08, right=0.85)
    ax[0,j].legend(bbox_to_anchor=(1.01, 0.5), fontsize=legend_font_size)

    output_filename = pines_path/('Objects/'+short_name+'/analysis/diagnostic_plots/'+short_name+'_cutout_positions.png')
    plt.savefig(output_filename, dpi=300)

    return

def absolute_image_position_plot(target, centroided_sources):
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)

    #Get plot style parameters. 
    title_size, axis_title_size, axis_ticks_font_size, legend_font_size = plot_style()
    
    #Get list of souce names in the centroid output.
    source_names = get_source_names(centroided_sources)
    centroided_sources.columns = centroided_sources.keys().str.strip()

    #Get times from the centroid output and split them by night. 
    times_full = np.array(centroided_sources['Time (JD UTC)'])
    night_inds = night_splitter(times_full)
    num_nights = len(night_inds)
    times_nights = [times_full[night_inds[i]] for i in range(num_nights)]
    standard_x = standard_x_range(times_nights)
    
    source = source_names[0]
    fig, ax = plt.subplots(nrows=2, ncols=num_nights, figsize=(17,9), sharex='col', sharey='row')
    for j in range(num_nights):
        if j == 0:
            ax[0,j].set_ylabel('Image X', fontsize=axis_title_size)
            ax[1,j].set_ylabel('Image Y', fontsize=axis_title_size)

        inds = night_inds[j]
        times = times_nights[j]
        absolute_x = np.array(centroided_sources[source+' Image X'][inds])
        absolute_y = np.array(centroided_sources[source+' Image Y'][inds])
        ax[0,j].plot(times, absolute_x, marker='.', linestyle='', alpha=0.3, color='tab:blue', label='Unbinned x position')
        ax[1,j].plot(times, absolute_y, marker='.', linestyle='', alpha=0.3, color='tab:orange', label='Unbinned y position')

        #bin
        block_inds = block_splitter(times)
        block_times= np.zeros(len(block_inds))
        block_x = np.zeros(len(block_inds))
        block_x_err = np.zeros(len(block_inds))
        block_y = np.zeros(len(block_inds))
        block_y_err = np.zeros(len(block_inds))
        for k in range(len(block_inds)):
            block_times[k] = np.mean(times[block_inds[k]])
            block_x[k] = np.mean(absolute_x[block_inds[k]])
            block_x_err[k] = np.std(absolute_x[block_inds[k]]) / np.sqrt(len(absolute_x[block_inds[k]]))
            block_y[k] = np.mean(absolute_y[block_inds[k]])
            block_y_err[k] = np.std(absolute_y[block_inds[k]]) / np.sqrt(len(absolute_y[block_inds[k]]))

        ax[0,j].errorbar(block_times, block_x, block_x_err, marker='o', linestyle='', color='tab:blue', ms=8, mfc='none', mew=2, label='Binned x position')
        ax[1,j].errorbar(block_times, block_y, block_y_err, marker='o', linestyle='', color='tab:orange', ms=8, mfc='none', mew=2, label='Binned y position')

        ax[0,j].tick_params(labelsize=axis_ticks_font_size)
        ax[1,j].tick_params(labelsize=axis_ticks_font_size)

        ax[0,j].grid(alpha=0.2)
        ax[1,j].grid(alpha=0.2)
        ax[1,j].set_xlabel('Time (JD UTC)', fontsize=axis_title_size)

    plt.suptitle(source+' Image Centroid Positions', fontsize=title_size)
    plt.subplots_adjust(left=0.07, hspace=0.05, wspace=0.05, top=0.92, bottom=0.08, right=0.85)
    ax[0,j].legend(bbox_to_anchor=(1.01, 0.5), fontsize=legend_font_size)
    ax[1,j].legend(bbox_to_anchor=(1.01, 0.5), fontsize=legend_font_size)
    ax[0,j].set_xlim(np.mean(times)-standard_x/2, np.mean(times)+standard_x/2)
    ax[1,j].set_xlim(np.mean(times)-standard_x/2, np.mean(times)+standard_x/2)

    output_filename = pines_path/('Objects/'+short_name+'/analysis/diagnostic_plots/'+source+'_image_positions.png')
    plt.savefig(output_filename, dpi=300)

    return 

def background_plot(target, centroided_sources, gain=8.21):
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)

    #Get plot style parameters. 
    title_size, axis_title_size, axis_ticks_font_size, legend_font_size = plot_style()

    analysis_path = pines_path/('Objects/'+short_name+'/analysis')
    phot_path = pines_path/('Objects/'+short_name+'/aper_phot')
    phot_files = np.array(natsorted([x for x in phot_path.glob('*.csv')]))
    if os.path.exists(analysis_path/('optimal_aperture.txt')):
        with open(analysis_path/('optimal_aperture.txt'), 'r') as f:
            best_ap = f.readlines()[0]
        ap_list = np.array([str(i).split('/')[-1].split('_')[3] for i in phot_files])
        best_ap_ind = np.where(ap_list == best_ap)[0][0]
    else:
        print('No optimal_aperture.txt file for {}.\nUsing first photometry file in {}.'.format(target, phot_path))
        best_ap_ind = 0

    phot_file = phot_files[best_ap_ind]
    phot_df = pines_log_reader(phot_file)

    backgrounds = np.array(phot_df[short_name+' Background'])/gain
    times_full = np.array(phot_df['Time JD'])
    night_inds = night_splitter(times_full)
    num_nights = len(night_inds)
    times_nights = [times_full[night_inds[i]] for i in range(num_nights)]
    standard_x = standard_x_range(times_nights)

    fig, ax = plt.subplots(nrows=1, ncols=num_nights, figsize=(17,5), sharey=True)
    for i in range(num_nights):
        if i == 0:
            ax[i].set_ylabel('Background (ADU)', fontsize=axis_title_size)

        inds = night_inds[i]
        ax[i].plot(times_full[inds], backgrounds[inds], marker='.', linestyle='', color='tab:orange', alpha=0.3, label='Unbinned background')
        ax[i].tick_params(labelsize=axis_ticks_font_size)
        ax[i].set_xlabel('Time (JD UTC)', fontsize=axis_title_size)
        ax[i].grid(alpha=0.2)
        ax[i].set_xlim(np.mean(times_full[inds])-standard_x/2, np.mean(times_full[inds]+standard_x/2))

        #bin
        block_inds = block_splitter(times_full[inds])
        block_x = np.zeros(len(block_inds))
        block_y = np.zeros(len(block_inds))
        block_y_err = np.zeros(len(block_inds))
        for j in range(len(block_inds)):
            block_x[j] = np.mean(times_full[inds][block_inds[j]])
            block_y[j] = np.mean(backgrounds[inds][block_inds[j]])
            block_y_err[j] = np.std(backgrounds[inds][block_inds[j]]) / np.sqrt(len(backgrounds[inds][block_inds[j]]))

        ax[i].errorbar(block_x, block_y, block_y_err, marker='o', linestyle='', color='tab:orange', ms=8, mfc='none', mew=2, label='Binned background')
        
        #Interpolate each night's seeing.
        fit_times = np.linspace(block_x[0], block_x[-1], 1000)
        interp = CubicSpline(block_x, block_y)
        interp_fit = interp(fit_times)
        ax[i].plot(fit_times, interp_fit, color='b', lw=2, zorder=0, alpha=0.7, label='CS Interpolation')

    ax[i].legend(bbox_to_anchor=(1.01, 0.5), fontsize=legend_font_size)
    plt.suptitle(short_name+' Background Measurements', fontsize=title_size)
    plt.subplots_adjust(left=0.07, hspace=0.05, wspace=0.05, top=0.92, bottom=0.17, right=0.85)

    output_filename = pines_path/('Objects/'+short_name+'/analysis/diagnostic_plots/'+short_name+'_backgrounds.png')
    plt.savefig(output_filename, dpi=300)
    return