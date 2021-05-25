import pdb
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check 
import pines_analysis_toolkit as pat 
from glob import glob 
from natsort import natsorted 
import numpy as np 

def align_and_stack_images(target): 
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    reduced_path = pines_path/('Objects/'+short_name+'/reduced')
    log_path = pines_path/'Logs/'

    reduced_files = np.array(natsorted([x for x in reduced_path.glob('*red.fits')]))
    dates = np.array(list(set([i.name.split('.')[0] for i in reduced_files])))
    log_names = np.array(natsorted([i+'_log.txt' for i in dates]))

    aligned_files = []
    for i in range(len(log_names)):
        log = log_path/log_names[i]
        log_df = pat.utils.pines_log_reader(log)
        locs = np.where((abs(log_df['X shift']) < 0.5) & (abs(log_df['Y shift']) < 0.5))[0]
        good_files = np.array(log_df['Filename'][np.where((abs(log_df['X shift']) < 0.5) & (abs(log_df['Y shift']) < 0.5))[0]])
        good_files = [i.split('.fits')[0]+'_red.fits' for i in good_files]
        aligned_files.extend(good_files)
    pdb.set_trace()
    return 

if __name__ == '__main__':
    target = '2MASS J05012406-0010452'
    align_and_stack_images(target)