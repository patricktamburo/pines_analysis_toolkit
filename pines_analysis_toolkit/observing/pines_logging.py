import numpy as np 
import pdb 

def pines_logging(filename, date, target_name, filter_name, exptime, airmass, x_shift, y_shift, x_seeing, y_seeing):
    '''
    Authors:
		Patrick Tamburo, Boston University, January 2021
	Purpose:
        Exports a line of log text in the same style as the telescope observing logs. 
	Inputs:
        filename (str): a targets 'long' 2MASS name, e.g. '2MASS J01234567+012345678'
        date (str): the UT date of the log whose shifts you want to update in YYYYMMDD format, e.g. '20151110'
        target_name (str): whether or not to push the updated log to the PINES server (only admins can do this)
        filter_name (str): filter name (J, H, K)
        exptime (str): exposure time in seconds
        airmass (str): airmass
        x_shift (str): x pixel shift
        y_shift (str): y pixel shift
        x_seeing (str): x seeing in arcsec
        y_seeing (str): y seeing in arcsec
    Outputs:
        log_text (str): a line of log text with the same spacing as those in the PINES telescope logs.
	TODO:
    FIXME:
    '''

    try:
        log_text = ' {:<19}, {:<20}, {:<30}, {:<6}, {:<8}, {:<8}, {:<8}, {:<8}, {:<9}, {:<9}\n'.format(filename, date, target_name, 
                                                                                                        filter_name,str(exptime),
                                                                                                        str(airmass),str(x_shift),
                                                                                                        str(y_shift),
                                                                                                        str(x_seeing),
                                                                                                        str(y_seeing))
    except:
        pdb.set_trace()

    return log_text