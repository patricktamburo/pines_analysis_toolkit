import os
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.utils.object_directory_creator import object_directory_creator
from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader
from pines_analysis_toolkit.data.get_master_log import get_master_log
from pines_analysis_toolkit.data.get_calibrations import get_calibrations
import numpy as np
import time
from progressbar import ProgressBar

def get_raw_science_files(sftp, target_name):
    """Downloads raw science files for a target

    :param sftp: sftp connection to the PINES server
    :type sftp: pysftp connection
    :param target_name: long name for the target
    :type target_name: str
    """

    t1 = time.time()
    print('')
    print('Starting get_raw_science files for {}.'.format(target_name))
   
    #Get the user's pines_analysis_toolkit path 
    pines_path = pines_dir_check()

    #Get the target's short name and set up a data directory, if necessary. 
    short_name = short_name_creator(target_name)
    if not os.path.exists(pines_path/('Objects/'+short_name)):
        object_directory_creator(pines_path, short_name)

    raw_data_path = pines_path/('Objects/'+short_name+'/raw/')
    dark_path = pines_path/('Calibrations/Darks')
    flats_path = pines_path/('Calibrations/Flats/Domeflats')
        
    #Grab an up-to-date copy of the master log, which will be used to find images. 
    get_master_log(sftp, pines_path)

    #Let's grab all of the available calibration data on pines.bu.edu.
    get_calibrations(sftp, pines_path)
    print('Calibrations up to date!')
    print('')
    time.sleep(2)

    #Read in the master target list and find images of the requested target. 
    df = pines_log_reader(pines_path/('Logs/master_log.txt'))
    targ_inds = np.where(np.array(df['Target']) == target_name)
    file_names = np.array(df['Filename'])[targ_inds]
    print('Searching pines.bu.edu for raw science files for {}.'.format(target_name))
    #Get list of dates that data are from, in chronological order. 
    dates = [int(i) for i in list(set([str.split(file_names[i],'.')[0] for i in range(len(file_names))]))]
    dates = np.array(dates)[np.argsort(dates)]
    comp_dates = np.array([int(file_names[i].split('.')[0]) for i in range(len(file_names))])
    print('Found ',len(file_names),' raw files for ',target_name,' on ',len(dates),' dates.')
    date_holder = [[] for x in range(len(dates))]
    for i in range(len(dates)):
        date = dates[i]
        print(date,': ',len(np.where(comp_dates == date)[0]),' files.')
        date_holder[i].extend(file_names[np.where(comp_dates == date)[0]])
        time.sleep(0.5)

    dates = [str(i) for i in dates]
    #Now download the identified data. 
    sftp.chdir('/data/raw/mimir')
    run_dirs = sftp.listdir()
    file_num = 1
    for i in range(len(run_dirs)):
        sftp.chdir(run_dirs[i])
        night_dirs = sftp.listdir()
        for j in range(len(night_dirs)):
            night_check = night_dirs[j]
            if night_check in dates:
                sftp.chdir(night_check)

                #Check that this night's log is in the /data/logs/ folder. If not, move a copy there.
                log_name = night_check+'_log.txt'
                if not log_name in sftp.listdir('/data/logs/'): 
                    sftp.get(log_name, pines_path/('Logs/'+log_name))
                    sftp.put(pines_path/('Logs/'+log_name), '/data/logs/'+log_name)
                
                date_holder_ind = np.where(np.array(dates) == night_check)[0][0]
                files = date_holder[date_holder_ind]
                print('Downloading {} data for {} to {}!'.format(night_check, target_name, raw_data_path))
                pbar = ProgressBar()
                for k in pbar(range(len(files))):
                    if not (raw_data_path/files[k]).exists():
                        #print('Downloading to {}, {} of {}'.format(raw_data_path/files[k],file_num,len(file_names)))
                        sftp.get(files[k],raw_data_path/files[k])
                    #else:
                        #print('{} already in {}, skipping download.'.format(files[k],raw_data_path))
                    file_num += 1
                sftp.chdir('..')
        sftp.chdir('..')
    
    print('')
    #Now grab the logs from /data/logs/ on the PINES server.
    sftp.chdir('/data/logs')
    for i in range(len(dates)):
        log_name = dates[i]+'_log.txt'
        
        print('Downloading {} to {}.'.format(log_name, pines_path/('Logs/'+log_name)))
        sftp.get(log_name, pines_path/('Logs/'+log_name))

        
    # sftp.chdir('/data/raw/mimir')
    # print('')
    # for i in range(len(run_dirs)):
    #         sftp.chdir(run_dirs[i])
    #         night_dirs = sftp.listdir()
    #         for j in range(len(night_dirs)):
    #             night_check = night_dirs[j]
    #             if night_check in dates:
    #                 sftp.chdir(night_check)
    #                 log_name = night_check+'_log.txt'
    #                 files_in_path = sftp.listdir()
    #                 if log_name in files_in_path:
    #                     if not (pines_path/('Logs/'+log_name)).exists():
    #                         sftp.get(log_name,pines_path/('Logs/'+log_name))
    #                         print('Downloading {} to {}.'.format(log_name, pines_path/('Logs/'+log_name)))
    #                     else:
    #                         print('{} already in {}, skipping download.'.format(log_name,pines_path/'Logs/'))
    #                 sftp.chdir('..')
    #         sftp.chdir('..')