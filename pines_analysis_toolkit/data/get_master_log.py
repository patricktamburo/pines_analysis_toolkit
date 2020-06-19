
import os
import pdb

'''Authors: 
        Patrick Tamburo, Boston University, June 2020
   Purpose: 
        Updates master_log.txt on the pines server, and downloads it if it doesn't match your local copy.
   Inputs:
        ssh (paramiko.client.SSHClient): ssh connection to pines.bu.edu.
        sftp (paramiko.sftp_client.SFTPClient): sftp connection to pines.pu.edu.
        pines_path (str): Path to local pines directory.
    Outputs:
        None
    TODO:
        Does this fail if you don't have write permissions? 
'''

def get_master_log(sftp, pines_path):
    #if (pines_path/('Logs/master_log.txt')).exists():

        #Run a remote command to update the master log. 
        #stdin, stdout, stderr = ssh.exec_command('/usr/local/bin/python3 /volume1/code/PINES_server_scripts/master_log_creator.py')

        #Count the lines in master_log.txt on the pines server, and only donwload it if it doesn't match your local master_log.txt
        #stdin, stdout, stderr = ssh.exec_command('wc -l /volume1/code/PINES_server_scripts/master_log.txt')
        #output = stdout.read()
        #remote_master_lines = int(str(output.splitlines()[0]).split(' ')[0].split("'")[1])

        #local_master_lines = int(len(open(pines_path/('Logs/master_log.txt')).readlines(  )))

        #f local_master_lines != remote_master_lines:
        #print('Downloading up-to-date copy of master_log.txt to {}/Logs/'.format(pines_path))
        #Download the master log
        #sftp.chdir('code/PINES_server_scripts/')
        #sftp.get('master_log.txt',(pines_path/'Logs/master_log.txt'))
        #else:
        #    print('master_log.txt up to date!')
    #else:

    print('Downloading master_log.txt to {}/Logs/'.format(pines_path))
    sftp.chdir('/code/PINES_server_scripts/')
    sftp.get('master_log.txt',(pines_path/'Logs/master_log.txt'))
    
    print('')
    return