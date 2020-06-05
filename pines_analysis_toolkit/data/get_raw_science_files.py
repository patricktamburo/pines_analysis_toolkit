import paramiko
import os
import pathlib
import pdb
from pines_analysis_pipeline.utils.pines_dir_check import pines_dir_check

'''Grabs raw science files from the PINES server.'''
def get_raw_science_files(username,password,run,date,local_path):
    #Check if Documents/PINES folder exists.
    pines_dir_check()
    
    #Open ssh connection and set up local/remote paths.
    ssh = paramiko.SSHClient()
    ssh.load_system_host_keys()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect('pines.bu.edu',username=username,password=password)
    a_path = pathlib.Path('/volume1/data/raw/mimir/'+run+'/'+date+'/')

    #Grab list of .fits files in the specified directory. 
    a_pattern = '*.fits'
    raw_command = 'find {path} -name {pattern}'
    command = raw_command.format(path=a_path,pattern=a_pattern)
    stdin, stdout, stderr = ssh.exec_command(command)
    file_list = stdout.read().splitlines()

    #Download the files. 
    sftp = ssh.open_sftp()
    for file in file_list: 
        (head,file_name) = os.path.split(file)
        print('Downloading ',file_name)
        pdb.set_trace()
        sftp.get(file, local_path)
        pdb.set_trace()
