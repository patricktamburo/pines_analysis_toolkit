from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
import numpy as np 
import matplotlib 
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt 
plt.ion()

def bad_shift_identifier(target, date, bad_shift_threshold=200.):
    """Identifies bad shift values in PINES observing logs. 

    :param target: name of the target as it is in the log
    :type target: str
    :param date: date of log to check (YYYYMMDD)
    :type date: str
    :param bad_shift_threshold: pixel threshold above which shifts are considered bad, defaults to 200
    :type bad_shift_threshold: float, optional
    """
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