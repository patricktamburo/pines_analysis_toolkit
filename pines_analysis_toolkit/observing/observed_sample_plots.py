import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pdb
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.pines_login import pines_login
from matplotlib.patches import Polygon
import os

def observed_sample_plots(upload=True):

    def plot_mwd(RA,Dec,observed_flag,org=0,title='Mollweide projection', projection='mollweide',observed_plot=0):
        ''' 
        Plots targets on the sky in a 'Mollweide' projection.
        RA, Dec are arrays of the same length.
        RA takes values in [0,360), Dec in [-90,90],
        which represent angles in degrees.
        org is the origin of the plot, 0 or a multiple of 30 degrees in [0,360).
        title is the title of the figure.
        projection is the kind of projection: 'mollweide', 'aitoff', 'hammer', 'lambert'
        '''

        x = np.remainder(RA+360-org,360) # shift RA values
        ind = x>180
        x[ind] -= 360    # scale conversion to [-180, 180]
        x=-x    # reverse the scale: East to the left
        x_tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210]) #Label in degrees
        #x_tick_labels = np.array([150,140,130,120,110,100,90,80,70,60,50,40,30,20,10,0,350,340,330,320,310,300,290,280,270,260,250,240,230,220,210]) #FinerLabel in degrees

        x_tick_labels = np.remainder(x_tick_labels+360+org,360)
        # x_tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])/15 #Label in hours
        # x_tick_labels = np.remainder(x_tick_labels+24+org/15,24)
        x_tick_labels = [int(i) for i in x_tick_labels]
        fig = plt.figure(figsize=(15*.8, 7*.8))
        ax = fig.add_subplot(111, projection=projection)
        #ax.scatter(np.radians(x),np.radians(Dec),color=color,alpha=0.4,zorder=1, label='Targets')  # convert degrees to radians
        for i in range(len(x)):
            if np.array(observed_flag)[i] == 0:
                color = 'k'
            else:
                color = 'k'
                if observed_plot == 1:
                    color = 'g' #Turn on observed targest plotting.
            ax.scatter(np.radians(x[i]),np.radians(Dec[i]),color=color,alpha=0.4,zorder=1,s=25)
        ax.set_yticklabels([str(int(i))+'$^\circ$' for i in np.round(ax.get_yticks()*180/np.pi)],fontsize=15)
        ax.title.set_fontsize(20)
        ax.set_xlabel('RA')
        ax.xaxis.label.set_fontsize(20)
        ax.set_ylabel("Dec")
        ax.yaxis.label.set_fontsize(20)
        ax.set_xticklabels([],fontsize=16)     # we add the scale on the x axis
        ax.grid(True,alpha=0.3)
        month_texts = ['Sep','Aug','Jul','Jun','May','Apr','Mar','Feb','Jan','Dec','Nov','Oct']
        for i in range(len(month_texts)):
            ax.text(-180*np.pi/180+15*np.pi/180+30*np.pi/180*i,-35*np.pi/180,month_texts[i],ha='center',va='center',fontsize=14)
        for i in range(len(x_tick_labels)):
            ax.text(-150*np.pi/180+30*np.pi/180*i,-22.5*np.pi/180,str(x_tick_labels[i])+'$^\circ$',ha='center',va='center',fontsize=15)
        
        #Plot monsoon season. 
        monsoon_x_vertices = np.array([-150,-150,-90,-90,-150])*np.pi/180
        monsoon_y_vertices = np.array([-90,90,90,-90,-90])*np.pi/180
        monsoon_polygon = Polygon(np.array([[monsoon_x_vertices[i],monsoon_y_vertices[i]] for i in range(len(monsoon_x_vertices))]),color='r',alpha=0.15,label='Flagstaff monsoon season')    
        ax.add_patch(monsoon_polygon)
        plt.show()
        return ax

    '''Plots the current sample as given in 'PINES sample.xlsx' on Google drive and uploads to the PINES website.'''
    pines_path = pines_dir_check()
    sample_path = pines_path/('Misc/PINES Sample.xlsx')
    print('Make sure an up-to-date copy of PINES Sample.xlsx exists in {}.'.format(pines_path/'Misc/'))
    print('Download from the PINES Google Drive.\n')

    df = pd.read_excel(sample_path)
    df = df.dropna(how='all') #Remove rows that are all NaNs. 

    good_locs = np.where(df['Good'] == 1)[0] #Get only "good" targets
    ra = np.array(df['RA (deg)'][good_locs])
    dec = np.array(df['Dec (deg)'][good_locs])
    group_ids = df['Group ID'][good_locs]
    observed_flag = df['Observed?'][good_locs]
    observed_groups = np.unique(np.array(group_ids)[np.where(observed_flag != 0)[0]]) #Get the groups that have been observed. 
    number_observed = len(np.array(group_ids)[np.where(observed_flag != 0)[0]])
    
    #Plot 1: Sky map of good targets based on group.
    print('Updating sky plot...')
    ax = plot_mwd(ra,dec,observed_flag,org=180,projection='mollweide',observed_plot=1)
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc=1, bbox_to_anchor=(1.1, 1.1), fontsize=16)
    ax.grid(alpha=0.2)

    group_id_inds = np.arange(0,max(group_ids)+1)

    #Now loop over group_id inds, and draw boundaries around each group.
    for i in group_id_inds:
        targs_in_group = np.where(group_ids == i)[0]
        try:
            cluster_coords = np.array([[ra[i],dec[i]] for i in targs_in_group])
        except:
            pdb.set_trace()
        hull = ConvexHull(cluster_coords)
        for s in range(len(hull.simplices)):
            simplex = hull.simplices[s]
            x = np.remainder(cluster_coords[simplex, 0]+360-180,360) # shift RA values
            ind = x>180
            x[ind] -= 360    # scale conversion to [-180, 180]
            x=-x    # reverse the scale: East to the left
            if i in observed_groups:
                color =  'g'
                ax.plot(x*np.pi/180, cluster_coords[simplex, 1]*np.pi/180,color = color,lw=2,zorder=0,alpha=0.6,label='Observed')
            else:
                color = 'k'
                ax.plot(x*np.pi/180, cluster_coords[simplex, 1]*np.pi/180,color = color,lw=2,zorder=0,alpha=0.6,label='Not yet observed')

    ax.grid(alpha=0.4)
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))

    ax.legend(by_label.values(), by_label.keys(), loc=1, bbox_to_anchor=(0.65, 0.225))
    ax.set_title('PINES sample \n '+str(int(max(group_ids)+1))+' groups, '+str(len(good_locs))+' targets',fontsize=20)
    plt.tight_layout()
    sky_map_output_path = pines_path/('Misc/updated_sky_plot.png')
    plt.savefig(sky_map_output_path,dpi=300)
    plt.close() 


    ntargs = len(df)
    #Now do magnitude/SpT histograms
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
        except: #Some values don't follow the normal +/- convention (they were upper limits in the Gagne sheet), so have to read them in differently. 
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
    fig,ax = plt.subplots(nrows=2,ncols=1,figsize=(18*scale_factor,15*scale_factor))
    bins = np.array([11.25,11.75,12.25,12.75,13.25,13.75,14.25,14.75,15.25,15.75, 16.25, 16.75]) - 0.25
    ax[0].hist(mags, bins=bins, histtype='step',lw=3, ls='--', label='Full sample')
    ax[0].hist(observed_mags, bins=bins, histtype='bar', label='Observed sample', color='tab:blue')
    ax[0].axvline(median_mag, color='r', label='Median $m_J$ = {:2.1f}'.format(median_mag))
    ticks = [11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5]
    ax[0].plot()
    ax[0].set_xticks(ticks)
    ax[0].set_xticklabels([str(i) for i in ticks])
    ax[0].set_xlabel('$m_J$',fontsize=20)
    ax[0].set_ylabel('Number of targets',fontsize=20)
    ax[0].tick_params(axis='both', which='major', labelsize=16)
    ax[0].legend(fontsize=16, loc='upper left')
    #ax[0].grid(alpha=0.2)

    ax[1].hist(SpT_number,bins=np.arange(-0.5,max(SpT_number)+0.5,1),histtype='step',lw=3,color='orange', ls='--', label='Full sample')
    ax[1].hist(observed_SpT_numbers,bins=np.arange(-0.5,max(SpT_number)+0.5,1),histtype='bar',lw=3,color='orange', label='Observed sample')
    ticks = np.arange(0,max(SpT_number),1)
    ax[1].set_xticks(ticks)
    ax[1].set_xticklabels(['L0','L1','L2','L3','L4','L5','L6','L7','L8','L9','T0','T1','T2','T3','T4','T5','T6','T7'])
    ax[1].set_xlabel('Spectral Type',fontsize=20)
    ax[1].set_ylabel('Number of targets',fontsize=20)
    ax[1].tick_params(axis='both', which='major', labelsize=16)
    ax[1].legend(fontsize=16, loc='upper right')
    #ax[1].grid(alpha=0.2)

    plt.tight_layout()
    histogram_output_path = pines_path/'Misc/target_histograms.png'
    plt.savefig(histogram_output_path, dpi=300)
    plt.close() 

    #Edit the observing.html page to update the number of observed targets. 
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
    
    edit_line_ind = np.where(['To date, PINES has observed' in i for i in lines])[0][0]
    edit_line = lines[edit_line_ind]
    edit_line = edit_line.replace(edit_line.split('<u>')[1].split('</u>')[0], str(number_observed))
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