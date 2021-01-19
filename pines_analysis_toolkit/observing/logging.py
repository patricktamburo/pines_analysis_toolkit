import numpy as np 

def logging(filename, date, target_name, filter_name, exptime, airmass, x_shift, y_shift, x_seeing, y_seeing):

    log_text = ' {:<19}, {:<20}, {:<30}, {:<6}, {:<8}, {:<8}, {:<8}, {:<8}, {:<9}, {:<9}\n'.format(filename, date, target_name, 
                                                                                                       filter_name,str(exptime),
                                                                                                       str(airmass),str(np.round(x_shift,1)),
                                                                                                       str(np.round(y_shift,1)),
                                                                                                       str(np.round(x_seeing,1)),
                                                                                                       str(np.round(y_seeing,1)))
    return log_text