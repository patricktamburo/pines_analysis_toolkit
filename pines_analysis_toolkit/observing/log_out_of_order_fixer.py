from pathlib import Path
from pines_analysis_toolkit.utils.pines_log_reader import pines_log_reader
from pines_analysis_toolkit.observing.pines_logging import pines_logging
import pdb 
import numpy as np
import os
from astropy.io import fits

def get_download_path():
    """Returns the default downloads path for linux or windows"""
    if os.name == 'nt':
        import winreg
        sub_key = r'SOFTWARE\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders'
        downloads_guid = '{374DE290-123F-4565-9164-39C4925E467B}'
        with winreg.OpenKey(winreg.HKEY_CURRENT_USER, sub_key) as key:
            location = winreg.QueryValueEx(key, downloads_guid)[0]
        return location
    else:
        return os.path.join(os.path.expanduser('~'), 'downloads')

def log_out_of_order_fixer(log_path, sftp):
    '''
    Authors:
		Patrick Tamburo, Boston University, January 2021
	Purpose:
        Fixes logs with out-of-order/duplicate filenames. This is a weird bug that can sometimes happen in our observing scripts.
	Inputs:
        log_path (pathlib.Path object): Path to the log. 
        sftp (pysftp connection): Connection to the PINES server. 
    Outputs:
        Writes corrected log file to disk.
	TODO:
    '''

    myfile = open(log_path, 'r')
    original_lines = myfile.readlines()
    myfile.close()

    #Remove any 'test.fits' lines, if they exist. 
    lines = []
    for i in range(len(original_lines)):
        if 'test.fits' in original_lines[i].split(',')[0]:
            continue
        else:
            lines.append(original_lines[i])
    
    #Sort everything based on file number. 
    log_filenums = []
    for i in range(len(lines)):
        if lines[i].strip()[0] == '#':
            log_filenums.append(0)
        else:
            log_filenums.append(int(lines[i].split(',')[0].split('.')[1]))

    sort_inds = np.argsort(log_filenums)
    lines = list(np.array(lines)[sort_inds])
    log_filenums = np.array(log_filenums)[sort_inds]

    num_del = 0
    #Check for lack of entries, multiple entries for this file number in the log. 

    for i in range(0, np.max(log_filenums)):
        file_num = i + 1
        entry_inds = np.where(log_filenums == file_num)[0]
        num_entries = len(entry_inds)

        #Check if there are no entries for this file number...if not, create one by finding the file on the PINES server,
        #downloading it, and reading its header. 
        if num_entries == 0:
            found = False

            missing_filename = log_path.name.split('_')[0]+'.'+str(file_num)+'.fits'
            if missing_filename != '20200206.607.fits': #corrupt file
                print('{} missing from log.'.format(missing_filename))
                print('Searching for file on PINES server...')

                pines_raw_path = '/data/raw/mimir/'
                runs = sftp.listdir(pines_raw_path)

                #The file will be in one of these directories if it's on the server
                run_guess = missing_filename[0:6]
                ind = np.where(np.array(runs) == run_guess)[0][0]

                if ind + 1 == len(runs):
                    inds = np.arange(ind-1, ind+1)
                    runs = np.array(runs)[inds]
                    runs = [runs[1], runs[0]]
                else:
                    inds = np.arange(ind-1, ind+2)
                    runs = np.array(runs)[inds]
                    runs = [runs[1], runs[0], runs[2]] 

                for jj in range(len(runs)):
                    nights = sftp.listdir(pines_raw_path+'/'+runs[jj])
                    nights = [i for i in nights if i[0] == '2']
                    for kk in range(len(nights)):
                        server_files = sftp.listdir(pines_raw_path+runs[jj]+'/'+nights[kk])
                        if missing_filename in server_files:
                            print('File found: {}'.format(pines_raw_path+runs[jj]+'/'+nights[kk]+'/'+missing_filename))
                            found = True
                            #Download the file and grab relevant parameters from the header. 
                            user_download_path = get_download_path()
                            sftp.get(pines_raw_path+runs[jj]+'/'+nights[kk]+'/'+missing_filename, user_download_path+'/'+missing_filename)
                            header = fits.open(user_download_path+'/'+missing_filename)[0].header
                            date = header['DATE']
                            if header['OBJECT'] == 'dummy':
                                print('ERROR: PINES_watchdog logged "dummy" for object filename; inspect the field and update its name manually.')
                            else:
                                target_name = header['OBJECT'].split('J')[0] +' J'+ header['OBJECT'].split('J')[1]
                            filter_name = header['FILTNME2']
                            exptime = str(header['EXPTIME'])
                            airmass = str(header['AIRMASS'])
                            #Use these guesses for shifts/seeings. Will be updated at a later step.
                            x_shift = str(0.0)
                            y_shift = str(0.0)
                            x_seeing = str(2.5)
                            y_seeing = str(2.5)

                            #Add the line to the log.
                            line = pines_logging(missing_filename, date, target_name, filter_name, exptime, airmass, x_shift, y_shift, x_seeing, y_seeing)
                            lines.insert(i+1, line)

                            #Delete the file from the Downloads folder. 
                            os.remove(user_download_path+'/'+missing_filename)
                            break
                
            if not found:
                print('{} not found in any raw directories on the PINES server... inspect manually.'.format(missing_filename))
                pdb.set_trace()

        #Check for multiple entries in the log...if so, delete the extras. 
        elif num_entries > 1:
            #Only retain the FIRST entry in the log. Assume later entries are mistakes. 
            num_del += 1
            for j in range(len(entry_inds)-1, 0, -1):
                del lines[entry_inds[j]]
                
        #Update the log on disk.
        with open(log_path, 'w') as f:
            for line in lines:
                f.write(line)

    return

if __name__ == '__main__':
    log_path = Path('/Users/tamburo/Documents/PINES_analysis_toolkit/Logs/20201002_log.txt')
    log_out_of_order_fixer(log_path)