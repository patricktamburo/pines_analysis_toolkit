import matplotlib
matplotlib.use('TkAgg')
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

def seeing_plot(target, centroided_sources, bin_mins='', force_output_path=''):
    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    #Get plot style parameters. 
    title_size, axis_title_size, axis_ticks_font_size, legend_font_size = plot_style()
    
    #Get list of souce names in the centroid output.
    source_names = get_source_names(centroided_sources)
    centroided_sources.columns = centroided_sources.keys().str.strip()

    #Read in times and seeing values. 
    times_full = np.array(centroided_sources['Time BJD TDB'], dtype='float')
    seeing = np.array(centroided_sources['Seeing'], dtype='float')

    #Split up times by nights.
    night_inds = night_splitter(times_full)
    num_nights = len(night_inds)
    times_nights = [times_full[night_inds[i]] for i in range(num_nights)]
    standard_x = standard_x_range(times_nights)

    fig, ax = plt.subplots(nrows=1, ncols=num_nights, figsize=(17,5), sharey=True)
    for i in range(num_nights):
        if i == 0:
            if num_nights == 1:
                ax.set_ylabel('Seeing (")', fontsize=axis_title_size)
            else:
                ax[i].set_ylabel('Seeing (")', fontsize=axis_title_size)

        inds = night_inds[i]
        

        #bin
        block_inds = block_splitter(times_nights[i], bin_mins=bin_mins)
        block_x = np.zeros(len(block_inds))
        block_y = np.zeros(len(block_inds))
        block_y_err = np.zeros(len(block_inds))
        for j in range(len(block_inds)):
            block_x[j] = np.nanmean(times_nights[i][block_inds[j]])
            block_y[j] = np.nanmean(seeing[inds][block_inds[j]])
            block_y_err[j] = np.nanstd(seeing[inds][block_inds[j]]) / np.sqrt(len(seeing[inds][block_inds[j]]))
       
        
        #Interpolate each night's seeing.
        fit_times = np.linspace(block_x[0], block_x[-1], 1000)
        interp = CubicSpline(block_x, block_y)
        interp_fit = interp(fit_times)
        
        if num_nights == 1:
            ax.plot(times_nights[i], seeing[inds], marker='.', linestyle='', alpha=0.3, label='Raw seeing')
            ax.tick_params(labelsize=axis_ticks_font_size)
            ax.set_xlabel('Time (BJD$_{TDB}$)', fontsize=axis_title_size)
            ax.grid(alpha=0.2)
            ax.set_xlim(np.mean(times_nights[i])-standard_x/2, np.mean(times_nights[i])+standard_x/2)
            ax.errorbar(block_x, block_y, block_y_err, marker='o', linestyle='', color='tab:blue', ms=8, mfc='none', mew=2, label='Bin seeing')
            ax.plot(fit_times, interp_fit, color='r', lw=2, zorder=0, alpha=0.7, label='CS Interp.')
        else:
            ax[i].plot(times_nights[i], seeing[inds], marker='.', linestyle='', alpha=0.3, label='Raw seeing')
            ax[i].tick_params(labelsize=axis_ticks_font_size)
            ax[i].set_xlabel('Time (BJD$_{TDB}$)', fontsize=axis_title_size)
            ax[i].grid(alpha=0.2)
            ax[i].set_xlim(np.mean(times_nights[i])-standard_x/2, np.mean(times_nights[i])+standard_x/2)
            ax[i].errorbar(block_x, block_y, block_y_err, marker='o', linestyle='', color='tab:blue', ms=8, mfc='none', mew=2, label='Bin seeing')
            ax[i].plot(fit_times, interp_fit, color='r', lw=2, zorder=0, alpha=0.7, label='CS Interp.')
    if num_nights == 1:
        ax.legend(bbox_to_anchor=(1.01, 0.5), fontsize=legend_font_size)
    else:
        ax[i].legend(bbox_to_anchor=(1.01, 0.5), fontsize=legend_font_size)
    plt.suptitle(short_name+' Seeing Measurements', fontsize=title_size)
    plt.subplots_adjust(left=0.07, wspace=0.05, top=0.92, bottom=0.17)

    output_filename = pines_path/('Objects/'+short_name+'/analysis/diagnostic_plots/'+short_name+'_seeing.png')
    plt.savefig(output_filename, dpi=300)
    return

def relative_cutout_position_plot(target, centroided_sources, force_output_path=''):
    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()
    short_name = short_name_creator(target)

    #Get plot style parameters. 
    title_size, axis_title_size, axis_ticks_font_size, legend_font_size = plot_style()
    
    #Get list of souce names in the centroid output.
    source_names = get_source_names(centroided_sources)
    centroided_sources.columns = centroided_sources.keys().str.strip()

    #Get times from the centroid output and split them by night. 
    times_full = np.array(centroided_sources['Time BJD TDB'])
    night_inds = night_splitter(times_full)
    num_nights = len(night_inds)
    times_nights = [times_full[night_inds[i]] for i in range(num_nights)]
    standard_x = standard_x_range(times_nights)

    #Get the box size (I don't like that this is being determined by using the mean of the data...output it from centroider?)    
    box_w = int(np.round(2*np.nanmean(np.array(centroided_sources['Reference 1 Cutout X'], dtype='float')),0))
        
    fig, ax = plt.subplots(nrows=2, ncols=num_nights, figsize=(17,9), sharey=True)
    plt.subplots_adjust(left=0.07, hspace=0.05, wspace=0.05, top=0.92, bottom=0.17)
    markers = ['+', 'x', '*', 'X']
    for j in range(num_nights):
        inds = night_inds[j]
        if j == 0: 
            if num_nights == 1:
                ax[0].set_ylabel('Cutout X Position', fontsize=axis_title_size)
                ax[1].set_ylabel('Cutout Y Position', fontsize=axis_title_size)
                ax[0].set_ylim(0.25,0.75)
            else:
                ax[0,j].set_ylabel('Cutout X Position', fontsize=axis_title_size)
                ax[1,j].set_ylabel('Cutout Y Position', fontsize=axis_title_size)
                ax[0,j].set_ylim(0.25,0.75)
        for i in range(len(source_names)):
            if i == 1:
                continue
            cutout_x = np.array(centroided_sources[source_names[i]+' Cutout X'][inds], dtype='float')
            cutout_y = np.array(centroided_sources[source_names[i]+' Cutout Y'][inds], dtype='float')

            if i == 0:
                marker = 'o'
                label = 'Target'
            else:
                marker = markers[(i-1) % len(markers)]
                label = 'Ref. '+str(i)

            if num_nights == 1:
                ax[0].plot(times_nights[j], cutout_x, marker=marker, label=label, linestyle='')
                ax[1].plot(times_nights[j], cutout_y, marker=marker, linestyle='')
            else:
                ax[0,j].plot(times_nights[j], cutout_x, marker=marker, label=label, linestyle='')
                ax[1,j].plot(times_nights[j], cutout_y, marker=marker, linestyle='')

        if num_nights == 1:
            ax[0].tick_params(labelsize=axis_ticks_font_size)
            ax[0].set_xticklabels([])
            ax[0].axhline(box_w/2, zorder=0, color='r', label='Center pix.', lw=2)
            ax[0].set_xlim(np.mean(times_nights[j])-standard_x/2, np.mean(times_nights[j])+standard_x/2)
            ax[0].grid(alpha=0.2)
            ax[1].tick_params(labelsize=axis_ticks_font_size)
            ax[1].axhline(box_w/2, zorder=0, color='r', label='Center pix.', lw=2)
            ax[1].set_xlim(np.mean(times_nights[j])-standard_x/2, np.mean(times_nights[j])+standard_x/2)
            ax[1].set_xlabel('Time (BJD$_{TDB}$)', fontsize=axis_title_size)
            ax[1].grid(alpha=0.2)
        else:
            ax[0,j].tick_params(labelsize=axis_ticks_font_size)
            ax[0,j].set_xticklabels([])
            ax[0,j].axhline(box_w/2, zorder=0, color='r', label='Center pix.', lw=2)
            ax[0,j].set_xlim(np.mean(times_nights[j])-standard_x/2, np.mean(times_nights[j])+standard_x/2)
            ax[0,j].grid(alpha=0.2)
            ax[1,j].tick_params(labelsize=axis_ticks_font_size)
            ax[1,j].axhline(box_w/2, zorder=0, color='r', label='Center pix.', lw=2)
            ax[1,j].set_xlim(np.mean(times_nights[j])-standard_x/2, np.mean(times_nights[j])+standard_x/2)
            ax[1,j].set_xlabel('Time (BJD$_{TDB}$)', fontsize=axis_title_size)
            ax[1,j].grid(alpha=0.2)

        if j == num_nights - 1:
            if num_nights == 1:
                ax[0].legend(bbox_to_anchor=(1.01, 1.0), fontsize=legend_font_size)
            else:
                ax[0,j].legend(bbox_to_anchor=(1.01, 1.0), fontsize=legend_font_size)

    plt.suptitle(short_name+' Cutout Centroid Positions', fontsize=title_size)

    output_filename = pines_path/('Objects/'+short_name+'/analysis/diagnostic_plots/'+short_name+'_cutout_positions.png')
    plt.savefig(output_filename, dpi=300)

    return

def absolute_image_position_plot(target, centroided_sources, force_output_path='', bin_mins=''):
    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()
    short_name = short_name_creator(target)

    #Get plot style parameters. 
    title_size, axis_title_size, axis_ticks_font_size, legend_font_size = plot_style()
    
    #Get list of souce names in the centroid output.
    source_names = get_source_names(centroided_sources)
    centroided_sources.columns = centroided_sources.keys().str.strip()

    #Get times from the centroid output and split them by night. 
    times_full = np.array(centroided_sources['Time BJD TDB'])
    night_inds = night_splitter(times_full)
    num_nights = len(night_inds)
    times_nights = [times_full[night_inds[i]] for i in range(num_nights)]
    standard_x = standard_x_range(times_nights)
    
    targ_ind = np.where(['2MASS' in i for i in source_names])[0][0]
    targ_name = source_names[targ_ind]
    fig, ax = plt.subplots(nrows=2, ncols=num_nights, figsize=(17,9), sharex='col', sharey='row')
    plt.subplots_adjust(left=0.07, hspace=0.05, wspace=0.05, top=0.92, bottom=0.17)
    for j in range(num_nights):
        if j == 0:
            if num_nights == 1:
                ax[0].set_ylabel('Image X', fontsize=axis_title_size)
                ax[1].set_ylabel('Image Y', fontsize=axis_title_size)
            else:
                ax[0,j].set_ylabel('Image X', fontsize=axis_title_size)
                ax[1,j].set_ylabel('Image Y', fontsize=axis_title_size)

        inds = night_inds[j]
        times = times_nights[j]
        absolute_x = np.array(centroided_sources[targ_name+' Image X'][inds], dtype='float')
        absolute_y = np.array(centroided_sources[targ_name+' Image Y'][inds], dtype='float')
        if num_nights == 1:
            ax[0].plot(times, absolute_x, marker='.', linestyle='', alpha=0.3, color='tab:blue', label='Raw x')
            ax[1].plot(times, absolute_y, marker='.', linestyle='', alpha=0.3, color='tab:orange', label='Raw y')
        else:
            ax[0,j].plot(times, absolute_x, marker='.', linestyle='', alpha=0.3, color='tab:blue', label='Raw x')
            ax[1,j].plot(times, absolute_y, marker='.', linestyle='', alpha=0.3, color='tab:orange', label='Raw y')

        #bin
        block_inds = block_splitter(times, bin_mins=bin_mins)
        block_times= np.zeros(len(block_inds))
        block_x = np.zeros(len(block_inds))
        block_x_err = np.zeros(len(block_inds))
        block_y = np.zeros(len(block_inds))
        block_y_err = np.zeros(len(block_inds))
        for k in range(len(block_inds)):
            try:
                block_times[k] = np.nanmean(times[block_inds[k]])
                block_x[k] = np.nanmean(absolute_x[block_inds[k]])
                block_x_err[k] = np.nanstd(absolute_x[block_inds[k]]) / np.sqrt(len(absolute_x[block_inds[k]]))
                block_y[k] = np.nanmean(absolute_y[block_inds[k]])
                block_y_err[k] = np.nanstd(absolute_y[block_inds[k]]) / np.sqrt(len(absolute_y[block_inds[k]]))
            except:
                pdb.set_trace()


        if num_nights == 1:
            ax[0].errorbar(block_times, block_x, block_x_err, marker='o', linestyle='', color='tab:blue', ms=8, mfc='none', mew=2, label='Bin x')
            ax[1].errorbar(block_times, block_y, block_y_err, marker='o', linestyle='', color='tab:orange', ms=8, mfc='none', mew=2, label='Bin y')
            ax[0].tick_params(labelsize=axis_ticks_font_size)
            ax[1].tick_params(labelsize=axis_ticks_font_size)
            ax[0].grid(alpha=0.2)
            ax[1].grid(alpha=0.2)
            ax[1].set_xlabel('Time (BJD$_{TDB}$)', fontsize=axis_title_size)
        else:
            ax[0,j].errorbar(block_times, block_x, block_x_err, marker='o', linestyle='', color='tab:blue', ms=8, mfc='none', mew=2, label='Bin x')
            ax[1,j].errorbar(block_times, block_y, block_y_err, marker='o', linestyle='', color='tab:orange', ms=8, mfc='none', mew=2, label='Bin y')
            ax[0,j].tick_params(labelsize=axis_ticks_font_size)
            ax[1,j].tick_params(labelsize=axis_ticks_font_size)
            ax[0,j].grid(alpha=0.2)
            ax[1,j].grid(alpha=0.2)
            ax[1,j].set_xlabel('Time (BJD$_{TDB}$)', fontsize=axis_title_size)

        if j == num_nights - 1:
            if num_nights == 1:
                ax[0].legend(bbox_to_anchor=(1.29, 1), fontsize=legend_font_size)
                ax[1].legend(bbox_to_anchor=(1.29, 1), fontsize=legend_font_size)
            else:
                ax[0,j].legend(bbox_to_anchor=(1.29, 1), fontsize=legend_font_size)
                ax[1,j].legend(bbox_to_anchor=(1.29, 1), fontsize=legend_font_size)

    plt.suptitle(targ_name+' Image Centroid Positions', fontsize=title_size)
    #plt.subplots_adjust(left=0.07, hspace=0.05, wspace=0.05, top=0.92, bottom=0.08, right=0.85)

    
   #ax.legend(bbox_to_anchor=(1.01, 1), fontsize=14)

    if num_nights == 1:
        ax[0].set_xlim(np.mean(times)-standard_x/2, np.mean(times)+standard_x/2)
        ax[1].set_xlim(np.mean(times)-standard_x/2, np.mean(times)+standard_x/2)
    else:
        ax[0,j].set_xlim(np.mean(times)-standard_x/2, np.mean(times)+standard_x/2)
        ax[1,j].set_xlim(np.mean(times)-standard_x/2, np.mean(times)+standard_x/2)

    output_filename = pines_path/('Objects/'+short_name+'/analysis/diagnostic_plots/'+targ_name+'_image_positions.png')
    plt.savefig(output_filename, dpi=300)

    return 

def background_plot(target, centroided_sources, gain=8.21, bin_mins='', force_output_path=''):

    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()
    short_name = short_name_creator(target)

    #Get plot style parameters. 
    title_size, axis_title_size, axis_ticks_font_size, legend_font_size = plot_style()

    analysis_path = pines_path/('Objects/'+short_name+'/analysis')
    phot_path = pines_path/('Objects/'+short_name+'/aper_phot')
    phot_files = np.array(natsorted([x for x in phot_path.glob('*.csv')]))

    if os.path.exists(analysis_path/('optimal_aperture.txt')):
        with open(analysis_path/('optimal_aperture.txt'), 'r') as f:
            best_ap = f.readlines()[0].split(':  ')[1].split('_')[0]
        ap_list = np.array([str(i).split('/')[-1].split('_')[4] for i in phot_files])
        best_ap_ind = np.where(ap_list == best_ap)[0][0]
    else:
        print('No optimal_aperture.txt file for {}.\nUsing first photometry file in {}.'.format(target, phot_path))
        best_ap_ind = 0

    phot_file = phot_files[best_ap_ind]
    phot_df = pines_log_reader(phot_file)

    backgrounds = np.array(phot_df[short_name+' Background'], dtype='float')/gain
    times_full = np.array(phot_df['Time BJD TDB'], dtype='float')
    night_inds = night_splitter(times_full)
    num_nights = len(night_inds)
    times_nights = [times_full[night_inds[i]] for i in range(num_nights)]
    standard_x = standard_x_range(times_nights)

    fig, ax = plt.subplots(nrows=1, ncols=num_nights, figsize=(17,5), sharey=True)
    for i in range(num_nights):
        if i == 0:
            if num_nights == 1:
                ax.set_ylabel('Background (ADU)', fontsize=axis_title_size)
            else:
                ax[i].set_ylabel('Background (ADU)', fontsize=axis_title_size)

        inds = night_inds[i]
        #bin
        block_inds = block_splitter(times_full[inds], bin_mins=bin_mins)
        block_x = np.zeros(len(block_inds))
        block_y = np.zeros(len(block_inds))
        block_y_err = np.zeros(len(block_inds))
        for j in range(len(block_inds)):
            block_x[j] = np.nanmean(times_full[inds][block_inds[j]])
            block_y[j] = np.nanmean(backgrounds[inds][block_inds[j]])
            block_y_err[j] = np.nanstd(backgrounds[inds][block_inds[j]]) / np.sqrt(len(backgrounds[inds][block_inds[j]]))

        block_x = block_x[~np.isnan(block_y)]
        block_y_err = block_y_err[~np.isnan(block_y)]
        block_y = block_y[~np.isnan(block_y)]
        
        #Interpolate each night's seeing.
        fit_times = np.linspace(block_x[0], block_x[-1], 1000)
        try:
            interp = CubicSpline(block_x, block_y)
        except:
            pdb.set_trace()
        interp_fit = interp(fit_times)

        if num_nights == 1:
            ax.plot(times_full[inds], backgrounds[inds], marker='.', linestyle='', color='tab:orange', alpha=0.3, label='Raw bkg.')
            ax.tick_params(labelsize=axis_ticks_font_size)
            ax.set_xlabel('Time (BJD$_{TDB}$)', fontsize=axis_title_size)
            ax.grid(alpha=0.2)
            ax.set_xlim(np.mean(times_full[inds])-standard_x/2, np.mean(times_full[inds]+standard_x/2))
            ax.errorbar(block_x, block_y, block_y_err, marker='o', linestyle='', color='tab:orange', ms=8, mfc='none', mew=2, label='Bin bkg.')
            ax.legend(bbox_to_anchor=(1.01, 0.5), fontsize=legend_font_size)
            ax.plot(fit_times, interp_fit, color='b', lw=2, zorder=0, alpha=0.7, label='CS Interp.')
        else:
            ax[i].plot(times_full[inds], backgrounds[inds], marker='.', linestyle='', color='tab:orange', alpha=0.3, label='Raw bkg.')
            ax[i].tick_params(labelsize=axis_ticks_font_size)
            ax[i].set_xlabel('Time (BJD$_{TDB}$)', fontsize=axis_title_size)
            ax[i].grid(alpha=0.2)
            ax[i].set_xlim(np.mean(times_full[inds])-standard_x/2, np.mean(times_full[inds]+standard_x/2))
            ax[i].errorbar(block_x, block_y, block_y_err, marker='o', linestyle='', color='tab:orange', ms=8, mfc='none', mew=2, label='Bin bkg.')
            if i == num_nights - 1:
                ax[i].legend(bbox_to_anchor=(1.01, 0.5), fontsize=legend_font_size)
            ax[i].plot(fit_times, interp_fit, color='b', lw=2, zorder=0, alpha=0.7, label='CS Interp.')

    plt.suptitle(short_name+' Background Measurements', fontsize=title_size)
    plt.subplots_adjust(left=0.07, wspace=0.05, top=0.92, bottom=0.17)
    
    output_filename = pines_path/('Objects/'+short_name+'/analysis/diagnostic_plots/'+short_name+'_backgrounds.png')
    plt.savefig(output_filename, dpi=300)
    return

def airmass_plot(target, centroided_sources, bin_mins='', force_output_path=''):
    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()
    short_name = short_name_creator(target)

    #Get plot style parameters. 
    title_size, axis_title_size, axis_ticks_font_size, legend_font_size = plot_style()

    #Get list of souce names in the centroid output.
    centroided_sources.columns = centroided_sources.keys().str.strip()
    airmasses = np.array(centroided_sources['Airmass'])
    times_full = np.array(centroided_sources['Time BJD TDB'])
    night_inds = night_splitter(times_full)
    num_nights = len(night_inds)
    times_nights = [times_full[night_inds[i]] for i in range(num_nights)]
    standard_x = standard_x_range(times_nights)

    fig, ax = plt.subplots(nrows=1, ncols=num_nights, figsize=(17,5), sharey=True)
    for i in range(num_nights):
        if i == 0:
            if num_nights == 1:
                ax.set_ylabel('Airmass', fontsize=axis_title_size)
            else:
                ax[i].set_ylabel('Airmass', fontsize=axis_title_size)

        inds = night_inds[i]

        #bin
        block_inds = block_splitter(times_full[inds], bin_mins=bin_mins)
        block_x = np.zeros(len(block_inds))
        block_y = np.zeros(len(block_inds))
        block_y_err = np.zeros(len(block_inds))
        for j in range(len(block_inds)):
            block_x[j] = np.mean(times_full[inds][block_inds[j]])
            block_y[j] = np.mean(airmasses[inds][block_inds[j]])
            block_y_err[j] = np.std(airmasses[inds][block_inds[j]]) / np.sqrt(len(airmasses[inds][block_inds[j]]))

        
        
        #Interpolate each night's seeing.
        fit_times = np.linspace(block_x[0], block_x[-1], 1000)
        interp = CubicSpline(block_x, block_y)
        interp_fit = interp(fit_times)

        if num_nights == 1:
            ax.plot(times_full[inds], airmasses[inds], marker='.', linestyle='', color='m', alpha=0.3, label='Raw airmass')
            ax.tick_params(labelsize=axis_ticks_font_size)
            ax.set_xlabel('Time (BJD$_{TDB}$)', fontsize=axis_title_size)
            ax.grid(alpha=0.2)
            ax.set_xlim(np.mean(times_full[inds])-standard_x/2, np.mean(times_full[inds]+standard_x/2))
            ax.errorbar(block_x, block_y, block_y_err, marker='o', linestyle='', color='m', ms=8, mfc='none', mew=2, label='Bin airmass')
            ax.plot(fit_times, interp_fit, color='c', lw=2, zorder=0, alpha=0.7, label='CS Interp.')
        else:
            ax[i].plot(times_full[inds], airmasses[inds], marker='.', linestyle='', color='m', alpha=0.3, label='Raw airmass')
            ax[i].tick_params(labelsize=axis_ticks_font_size)
            ax[i].set_xlabel('Time (BJD$_{TDB}$)', fontsize=axis_title_size)
            ax[i].grid(alpha=0.2)
            ax[i].set_xlim(np.mean(times_full[inds])-standard_x/2, np.mean(times_full[inds]+standard_x/2))
            ax[i].errorbar(block_x, block_y, block_y_err, marker='o', linestyle='', color='m', ms=8, mfc='none', mew=2, label='Bin airmass')
            ax[i].plot(fit_times, interp_fit, color='c', lw=2, zorder=0, alpha=0.7, label='CS Interp.')
    
    if num_nights == 1:
        ax.legend(bbox_to_anchor=(1.005, 0.5), fontsize=legend_font_size)
    else:
        ax[i].legend(bbox_to_anchor=(1.005, 0.5), fontsize=legend_font_size)
    plt.suptitle(short_name+' Airmass Measurements', fontsize=title_size)
    plt.subplots_adjust(left=0.07, wspace=0.05, top=0.92, bottom=0.17)

    output_filename = pines_path/('Objects/'+short_name+'/analysis/diagnostic_plots/'+short_name+'_airmasses.png')
    plt.savefig(output_filename, dpi=300)
    return