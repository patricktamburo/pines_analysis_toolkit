import pdb
import os
import stat

'''Authors:
		Patrick Tamburo, Boston University, June 2020
	Purpose:
        Pulls all master dome flats from pines.bu.edu:/data/calibrations/Flats/Domeflats/
	Inputs:
		sftp (paramiko.sftp_client.SFTPClient): the sftp client to pines.bu.edu
        pines_dark_path (str): path to the PINES_analysis_toolkit/Calibrations/Darks/ directory
	Outputs:
		None
	TODO:
		None
'''
def get_master_darks(sftp, pines_dark_path):
    sftp.chdir('/data/calibrations/Darks')
    darks = sftp.listdir()
    for i in range(len(darks)):
        dark = darks[i]
        if not (pines_dark_path/('Master Darks/'+dark)).exists():
            sftp.get(dark, pines_dark_path/('Master Darks/'+dark))
            print('Downloading to {} '.format(pines_dark_path/('/Master Darks/'+dark)))
    sftp.chdir('..')
    sftp.chdir('..')
    sftp.chdir('..')
    sftp.chdir('..')
    return
