import pdb 

def pines_logging(filename, date, target_name, filter_name, exptime, airmass, x_shift, y_shift, x_seeing, y_seeing, post_processing_flag,   shift_quality_flag):
    """Exports a line of log text in the same style as the telescope observing logs. 

    :param filename: string of raw filename (e.g., '20210214.203.fits')
    :type filename: str
    :param date: date of image (YYYYMMDD)
    :type date: str
    :param target_name: taget's long name
    :type target_name: str
    :param filter_name: filter name (e.g., 'J' or 'H')
    :type filter_name: str
    :param exptime: exposure time
    :type exptime: str
    :param airmass: airmass
    :type airmass: str
    :param x_shift: x shift in pixels
    :type x_shift: str
    :param y_shift: y shift in pixels
    :type y_shift: str
    :param x_seeing: x seeing in arcsec
    :type x_seeing: str
    :param y_seeing: y seeing in arcsec
    :type y_seeing: str
    :param post_processing_flag: 0 (not yet post-processed) or 1 (post-processed)
    :type post_processing_flag: int
    :param shift_quality_flag: 0 (good shifts) or 1 (bad shifts)
    :type shift_quality_flag: int
    """
    try:
        log_text = ' {:<19}, {:<20}, {:<30}, {:<6}, {:<8}, {:<8}, {:<8}, {:<8}, {:<9}, {:<7}, {:<21}, {:<20}\n'.format(filename, date, target_name, 
                                                                                                        filter_name,str(exptime),
                                                                                                        str(airmass),str(x_shift),
                                                                                                        str(y_shift),
                                                                                                        str(x_seeing),
                                                                                                        str(y_seeing),
                                                                                                        str(post_processing_flag),
                                                                                                        str(shift_quality_flag))
    except:
        pdb.set_trace()

    return log_text