import pdb
import os
import stat

'''Authors:
		Patrick Tamburo, Boston University, June 2020
	Purpose:
        Pulls all master dome flats from pines.bu.edu:/data/calibrations/Flats/Domeflats/
	Inputs:
		sftp (paramiko.sftp_client.SFTPClient): the sftp client to pines.bu.edu
        pines_dome_flat_path (pathlib.Path object): path to the PINES_analysis_toolkit/Calibrations/Flats/Domeflats directory
    Outputs:
		None
	TODO:
		None
'''
def get_master_dome_flats(sftp, pines_dome_flat_path):
    sftp.chdir('/data/calibrations/Flats/Domeflats')
    bands = sftp.listdir()
    for i in range(len(bands)):
        band = bands[i]
        sftp.chdir(band)
        flats = sftp.listdir()
        for j in range(len(flats)):
            flat = flats[j]
            if not (pines_dome_flat_path/(band+'/Master Flats/'+flat)).exists():
                sftp.get(flat, pines_dome_flat_path/(band+'/Master Flats/'+flat))
                print('Downloading to {}'.format(pines_dome_flat_path/(band+'/Master Flats/'+flat)))
        sftp.chdir('..')
    sftp.chdir('..')
    sftp.chdir('..')
    sftp.chdir('..')
    sftp.chdir('..')
    return