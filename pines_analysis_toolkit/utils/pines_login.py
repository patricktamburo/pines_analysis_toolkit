import pysftp
import getpass
import pdb

'''Authors: 
        Patrick Tamburo, Boston University, June 2020
   Purpose: 
        Logs user into the pines.bu.edu server with an sftp connection. 
   Inputs:
        None
   Outputs:
        sftp (pysftp.Connection): the sftp connection to the pines server. 
    TODO:
        None 
'''

def pines_login():
    i = 3
    while i != 0:
        print('')
        username = input('Enter pines.bu.edu username: ')
        password = getpass.getpass('Enter password for {}@pines.bu.edu: '.format(username))
        print('')
        try:
            sftp = pysftp.Connection('pines.bu.edu', username=username, password=password)
            print('Login successful!\n')
            return sftp
        except:
            i -= 1
            if i == 1:
                verb = 'try' #lol
            else:
                verb = 'tries'
            if i != 0:
                print('Login failed! {} {} remaining.'.format(i, verb))
            else:
                print('ERROR: login failed after {} tries.'.format(i))