import pdb
import time
def short_name_creator(long_name):
    """Tries to create a 'short name' out of a full 2MASS identifier. 
    I.e., 2MASS J00001354+2554180 would be converted to 2MASS 0000+2554. 

    :param long_name: full 2MASS name of the target
    :type long_name: str
    :return: short name for the target 
    :rtype: str
    """
    name = long_name.replace(' ','')

    if '2MASS' in name:
        #If the user already passed a 2MASS short name (i.e., 2MASS XXXX+YYYY), just return it.
        if len(name) == 14:
            short_name = long_name
        
        #Otherwise, create the short name.
        else:
            if '+' in name:
                short_name = '2MASS '+name.split('J')[1].split('+')[0][0:4]+'+'+name.split('J')[1].split('+')[1][0:4]
            elif '-' in name:
                short_name = '2MASS '+name.split('J')[1].split('-')[0][0:4]+'-'+name.split('J')[1].split('-')[1][0:4]
    #If a name other than a 2MASS identifier was passed, just return the long name. 
    else:
        print('')
        print('WARNING: long_name does not match 2MASS Jxxxxxxxx+xxxxxxx format, returning long_name.')
        short_name = long_name
    
    return short_name