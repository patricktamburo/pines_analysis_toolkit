import paramiko
import os
import pathlib
import pdb
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.utils.object_directory_creator import object_directory_creator
from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader
from pines_analysis_toolkit.data.get_master_log import get_master_log
import getpass 
import pandas
import numpy as np
import time

'''Authors: 
        Patrick Tamburo, Boston University, June 2020
   Purpose: 
        Finds raw science files on the PINES server for a specified target, and downloads them.
   Inputs:
        target_name (str): the target's full 2MASS name, i.e. '2MASS J01234567+0123456' 
    Outputs:
        None
    TODO:
        Zip files for faster download
        Allow user to specify run if target is observed during multiple runs 
        Have program automatically regenerate master_log.txt on the PINES server when run
'''

def get_raw_science_files(target_name):
    #Prompt login: 
    print('')
    username = input('Enter username: ')
    password = getpass.getpass('Enter password: ')
    print('')


    #Get the user's pines_analysis_toolkit path 
    pines_path = pines_dir_check()

    #Get the target's short name and set up a data directory, if necessary. 
    short_name = short_name_creator(target_name)
    if not os.path.exists(pines_path/('Objects/'+short_name)):
        object_directory_creator(pines_path, short_name)

    raw_data_path = pines_path/('Objects/'+short_name+'/raw/')

    #Open ssh connection and set up local/remote paths.
    ssh = paramiko.SSHClient()
    ssh.load_system_host_keys()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect('pines.bu.edu',username=username, password=password)
    sftp = ssh.open_sftp()

    t1 = time.time()
    print('Searching pines.bu.edu for raw science files for {}.'.format(target_name))

    username = ''
    password = ''
    
    #Grab an up-to-date copy of the master log, which will be used to find images. 
    get_master_log(ssh, sftp, pines_path)

    #Read in the master target list and find images of the requested target. 
    df = pines_log_reader(pines_path/('Logs/master_log.txt'))
    targ_inds = np.where(np.array(df['Target']) == target_name)
    file_names = np.array(df['Filename'])[targ_inds]
    print('')
    
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
    print('')

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
                date_holder_ind = np.where(np.array(dates) == night_check)[0][0]
                files = date_holder[date_holder_ind]
                for k in range(len(files)):
                    print('Downloading to {}, {} of {}'.format(raw_data_path/files[k],file_num,len(file_names)))
                    if not (raw_data_path/files[k]).exists():
                        sftp.get(files[k],raw_data_path/files[k])
                    else:
                        print('{} already in {}, skipping.'.format(files[k],raw_data_path))
                    file_num += 1
                sftp.chdir('..')
        sftp.chdir('..')
        
    sftp.close()
    print('')
    print('get_raw_science_files runtime: ', np.round((time.time()-t1)/60,1), ' minutes.')
    print('Done!')
