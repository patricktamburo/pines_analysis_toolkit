import numpy as np 
import natsort 
import pines_analysis_toolkit as pat 
from progressbar import ProgressBar
import pdb 

def upload_reduced_data(sftp, short_name): 
    """Standalone function for uploading reduced data to the PINES server

    :param sftp: sftp connection to the PINES server
    :type sftp: pysftp connection
    :param short_name: short name for the target
    :type short_name: str 
    """
    pines_path = pat.utils.pines_dir_check()
    reduced_path = pines_path/('Objects/'+short_name+'/reduced')

    print('')
    print('Beginning upload process to pines.bu.edu...')
    print('NOTE:    Only PINES admins are able to upload.')
    print('WARNING: If these reduced images already exist on the PINES server, they will be overwritten!')

    sftp.chdir('/data/reduced/mimir/')
    files_to_upload = np.array(natsort.natsorted(np.array([x for x in reduced_path.glob('*red.fits')])))
    print('Uploading reduced {} data to the PINES server!'.format(short_name))
    pbar = ProgressBar()
    for i in pbar(range(len(files_to_upload))):
        file = files_to_upload[i]
        night_name = files_to_upload[i].name.split('.')[0]
        for dir_change in sftp.listdir():
            sftp.chdir(dir_change)
            nights = sftp.listdir()
            if night_name in nights:
                ind  = np.where(np.array(nights) == night_name)[0][0]
                break
            sftp.chdir('..')
        if nights[ind] != night_name:
            print('ERROR: the date of the file you want to upload does not match the date directory where the program wants to upload it.')
            pdb.set_trace()
        
        #print('Uploading to {}/{}, {} of {}'.format(sftp.getcwd(),nights[ind]+'/'+file.name, i+1,len(files_to_upload)))
        sftp.put(file,nights[ind]+'/'+file.name)
        sftp.chdir('..')
        