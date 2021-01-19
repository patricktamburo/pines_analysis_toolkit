
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
    print('Downloading master_log.txt to {}/Logs/'.format(pines_path))
    sftp.chdir('/code/PINES_server_scripts/')
    sftp.get('master_log.txt',(pines_path/'Logs/master_log.txt'))
    return