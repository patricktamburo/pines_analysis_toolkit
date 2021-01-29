from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
import pdb 
import numpy as np 
import matplotlib 
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt 
plt.ion()
from astropy.stats import sigma_clipped_stats

def bad_shift_identifier(target, date, bad_shift_threshold=200):
    pines_path = pines_dir_check()
    log_path = pines_path/('Logs/'+date+'_log.txt')
    log = pines_log_reader(log_path)
    target_inds = np.where(log['Target'] == target)[0]
    x_shifts = np.array(log['X shift'][target_inds])
    y_shifts = np.array(log['Y shift'][target_inds])

    bad_shift_inds = np.where((x_shifts > bad_shift_threshold) | (y_shifts > bad_shift_threshold))[0]
    shift_flags = np.zeros(len(target_inds), dtype=int)
    shift_flags[bad_shift_inds] = 1
    

    return 

if __name__ == '__main__':
    target = 'SIMP 0136'
    date = '20151110'
    bad_shift_identifier(target, date)