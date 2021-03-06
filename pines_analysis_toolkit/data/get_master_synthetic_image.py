import pysftp
import pdb 
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check

def get_master_synthetic_image(sftp, target_name):
    pines_path = pines_dir_check()
    sftp.chdir('/data/master_images/')
    synthetic_filename = target_name.replace(' ', '')+'_master_synthetic.fits'
    download_path = pines_path/('Calibrations/Master Synthetic Images/'+synthetic_filename)
    sftp.get(synthetic_filename, download_path)
    print('Downloaded {} to {}!'.format(synthetic_filename, download_path))