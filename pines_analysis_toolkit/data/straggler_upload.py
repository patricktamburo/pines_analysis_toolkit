import pdb 
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.utils.pines_login import pines_login
from natsort import natsorted 
from glob import glob 
import numpy as np 

def straggler_upload(sftp, targets, delete_raw=False):
    sftp.chdir('/data/reduced/mimir/')
    runs = sftp.listdir()
    pines_path = pines_dir_check()
    for i in range(len(targets)):
        target = targets[i]
        short_name = short_name_creator(target)
        raw_path = pines_path/('Objects/'+short_name+'/raw/')
        red_path = pines_path/('Objects/'+short_name+'/reduced/')
        raw_files = np.array(natsorted(glob(str(raw_path)+'/*.fits')))
        for j in range(len(raw_files)):
            sftp.chdir('/data/reduced/mimir/')
            red_file = str(red_path)+'/'+raw_files[j].split('/')[-1].split('.fits')[0]+'_red.fits'
            night = red_file.split('/')[-1].split('.')[0]
            run_guess = night[0:6]
            ind = np.where(np.array(runs) == run_guess)[0][0]
            if ind + 1 == len(runs):
                inds = np.arange(ind-1, ind+1)
                runs = np.array(runs)[inds]
                runs = [runs[1], runs[0]]
            else:
                inds = np.arange(ind-1, ind+2)
                runs = np.array(runs)[inds]
                runs = [runs[1], runs[0], runs[2]] 
            
            found = False
            for k in range(len(runs)):
                sftp.chdir(runs[k])
                if night in sftp.listdir():
                    print('Uploading {} to {}.'.format(red_file.split('/')[-1], sftp.pwd+'/'+night))
                    sftp.chdir(night)
                    sftp.put(red_file, red_file.split('/')[-1])
                    found = True
                else:
                    sftp.chdir('..')
                if found:
                    break
                pdb.set_trace()

if __name__ == '__main__':
    targets =  ['2MASS J11250502+0556424',
                '2MASS J11394192-0310039',
                '2MASS J11491231-0153003',
                '2MASS J11555389+0559577',
                '2MASS J11582077+0435014']
    sftp = pines_login()
    straggler_upload(sftp, targets)