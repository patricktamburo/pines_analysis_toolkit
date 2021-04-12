import pines_analysis_toolkit as pat 

def log_repair(sftp): 
    #Download all logs. 
    sftp.chdir('/data/logs')
    for i in range(len(dates)):
        log_name = dates[i]+'_log.txt'
        
        print('Downloading {} to {}.'.format(log_name, pines_path/('Logs/'+log_name)))
        sftp.get(log_name, pines_path/('Logs/'+log_name))

        pdb.set_trace()
    return 

if __name__ == '__main__':
    sftp = pat.utils.pines_login()
    log_repair(sftp)
