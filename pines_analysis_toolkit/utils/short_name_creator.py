import pdb
import time
def short_name_creator(long_name):
    name = long_name.replace(' ','')

    if '+' in name:
        short_name = '2MASS '+name.split('J')[1].split('+')[0][0:4]+'+'+name.split('J')[1].split('+')[1][0:4]
    elif '-' in name:
        short_name = '2MASS '+name.split('J')[1].split('-')[0][0:4]+'-'+name.split('J')[1].split('-')[1][0:4]
    else:
        print('')
        print('WARNING: long_name does not match 2MASS Jxxxxxxxx+xxxxxxx format, returning long_name.')
        #time.sleep(3)
        short_name = long_name
    
    return short_name