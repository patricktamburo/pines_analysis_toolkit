from astropy.io import fits 


def update_header(file_path, header_key, new_value):
    ''' 
        Author: 
            Patrick Tamburo, February 2021
        Purpose: 
            Updates a fits header key with a chosen value.
        Inputs:
            path (str): path to fits file
            header_key (str): the header key you want to update in this file (e.g., 'OBJECT', 'FILTNME2', etc.)
            new_value (flexible): the value you want to update the header key for this file to
        Outputs:
            Updates header of fits file.
        TODO:
    '''

    fits.setval(file_path, header_key, value=new_value)

