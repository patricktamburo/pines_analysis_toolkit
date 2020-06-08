import pdb
def short_name_creator(long_name):
    long_name = long_name.replace(' ','')
    if '+' in long_name:
        short_name = '2MASS '+long_name.split('J')[1].split('+')[0][0:4]+'+'+long_name.split('J')[1].split('+')[1][0:4]
    elif '-' in long_name:
        short_name = '2MASS '+long_name.split('J')[1].split('-')[0][0:4]+'-'+long_name.split('J')[1].split('-')[1][0:4]
    else:
        print('ERROR: Is this long_name correct?')
    return short_name