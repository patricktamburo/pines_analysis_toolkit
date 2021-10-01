from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader
from glob import glob
import pdb 
from natsort import natsorted 
import numpy as np 

def average_seeing(log_path):
    """Calculates average seeing from values in the PINES observing logs. 

    :param log_path: ~/PINES_analysis_toolkit/Logs/ path
    :type log_path: pathlib.PosixPath
    :return: average seeing 
    :rtype: float
    """
    try:
        df = pines_log_reader(log_path)
        if 'X seeing' in df.keys():
            seeing = np.array(df['X seeing'], dtype=float)
            seeing = np.array(seeing[np.where(~np.isnan(seeing))], dtype=float)
            seeing = seeing[np.where((seeing > 1.2) & (seeing < 7.0))[0]]
            mean_seeing = np.nanmean(seeing)
            std_seeing = np.nanstd(seeing)
            print('Average seeing for {}: {:1.1f} +/- {:1.1f}"'.format(log_path.split('/')[-1].split('_')[0], mean_seeing, std_seeing))
            return mean_seeing
    except:
        print('{}: No seeing measurements, inspect manually.'.format(log_path.split('/')[-1].split('_')[0]))
        return np.nan
    

if __name__ == '__main__':
    pines_dir = pines_dir_check()
    log_path = pines_dir/'Logs'
    logs = np.array(natsorted(glob(str(log_path/'*.txt'))))

    #Remove one bad log 
    bad_loc = np.where(np.array([logs[i].split('/')[-1] for i in range(len(logs))]) == '20201003_log.txt')[0][0]
    logs = np.delete(logs, bad_loc)
    #logs = logs[6:-1] #Ignore first set of logs, ignore master_log

    
    seeings = np.zeros(len(logs))
    for i in range(len(logs)):
        log_path = logs[i]
        seeings[i] = average_seeing(log_path)
    
    seeings = seeings[seeings < 5] #Cut some bad values. 
    print('Global average seeing: {:1.1f} +/- {:1.1f}.'.format(np.nanmean(seeings), np.nanstd(seeings)))
    pdb.set_trace()