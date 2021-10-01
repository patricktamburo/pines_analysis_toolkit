import pines_analysis_toolkit as pat 
import pdb 

def make_calibrations(sftp, exptimes, bands, dark_dates, flat_dates, dark_starts=[], dark_stops=[], lights_on_starts=[],
                      lights_on_stops=[], lights_off_starts=[], lights_off_stops=[], dark=True, flat=True, variable_pixels=True, 
                      hot_pixels=True, dead_pixels=True, bpm=True):
    """Makes all calibration files for a target. 

    :param sftp: sftp connection to the PINES server
    :type sftp: pysftp connection
    :param exptimes: list of exposure times for this target. E.g., if you observed with 15 and 30 seconds, pass [15, 30]
    :type exptimes: list of floats
    :param bands: list of bands used to observe this target. E.g., if you observed in J and H bands, pass ['J', 'H']. 
    :type bands: list of str
    :param dark_dates: list of dates on which darks were taken (order must match order of exptimes list, and there should be one entry for every entry in exptimes!)
    :type dark_dates: list of str
    :param flat_dates: list of dates on whic flats were taken (order must match order of bands list, and there should be one entry for every entry in exptimes!)
    :type flat_dates: list of str
    :param dark_starts: list of the start filenumbers for darks (order matches order of exptimes!), defaults to []
    :type dark_starts: list of int, optional
    :param dark_stops: list of the stop filenumbers for darks (order matches order of exptimes!), defaults to []
    :type dark_stops: list of int, optional
    :param lights_on_starts: list of the start filenumbers for lights on flats (order matches order of bands!), defaults to []
    :type lights_on_starts: list of int, optional
    :param lights_on_stops: list of the stop filenumbers for lights on flats (order matches order of bands!), defaults to []
    :type lights_on_stops: list of int, optional
    :param lights_off_starts: list of the start filenumbers for lights off flats (order matches order of bands!), defaults to []
    :type lights_off_starts: list of int, optional
    :param lights_off_stops: list of the stop filenumbers for lights off flats (order matches order of bands!), defaults to []
    :type lights_off_stops: list of int, optional
    :param dark: Whether or not to make the master darks, defaults to True
    :type dark: bool, optional
    :param flat: Whether or not to make the master flats, defaults to True
    :type flat: bool, optional
    :param variable_pixels: Whether or not to make the variable pixel masks, defaults to True
    :type variable_pixels: bool, optional
    :param hot_pixels: Whether or not to make the hot pixel masks, defaults to True
    :type hot_pixels: bool, optional
    :param dead_pixels: Whether or not to make the dead pixel masks, defaults to True
    :type dead_pixels: bool, optional
    :param bpm: Whether or not to make the bad pixel masks, defaults to True
    :type bpm: bool, optional
    """        

    pines_path = pat.utils.pines_dir_check()
    calibration_path = pines_path/('Calibrations/')

    assert len(bands) == len(flat_dates), 'Num. of elements in the bands list does not equal the num. of elements of the flat_dates list!\n You must have one flat_dates entry for every bands entry.'
    assert len(exptimes) == len(dark_dates), 'Num. of elements in the exptimes list does not equal the num. of elements of the dark_dates list!\n You must have one dark_dates entry for every exptimes entry.'

    #Master darks. 
    if dark:
        for i in range(len(exptimes)):
            dark_date = dark_dates[i]
            exptime = float(exptimes[i])
            dark_path = calibration_path/('Darks/Master Darks/master_dark_'+str(exptime)+'_s_'+dark_date+'.fits')
            #If the master dark for this exptime/date already exist, skip! 
            if dark_path.exists():
                print('{} already exists in {}, skipping.'.format(dark_path.name, calibration_path))
            else:
                print('Making master dark for {}-s exposure time on {}.'.format(exptime, dark_date))
                if len(dark_starts) == 0:
                    dark_start = 0
                    dark_stop  = 0
                else:
                    dark_start = dark_starts[i]
                    dark_stop  = dark_stops[i]
                pat.data.dark(dark_date, exptime, upload=True, delete_raw=True, sftp=sftp, dark_start=dark_start, dark_stop=dark_stop)

    print('')

    #Master flats. 
    if flat:
        for i in range(len(bands)):
            flat_date = flat_dates[i]
            band = bands[i]
            flat_path = calibration_path/('Flats/Domeflats/'+band+'/Master Flats/master_flat_'+band+'_'+flat_date+'.fits')
            if flat_path.exists():
                print('{} already exists in {}, skipping.'.format(flat_path.name, calibration_path))
            else:
                print('Making master flat for {}-band on {}.'.format(band, dark_date))

                if len(lights_on_starts) == 0:
                    lights_on_start  = 0
                    lights_on_stop   = 0
                    lights_off_start = 0
                    lights_off_stop  = 0
                else:
                    lights_on_start  = lights_on_starts[i]
                    lights_on_stop   = lights_on_stops[i]
                    lights_off_start = lights_off_starts[i]
                    lights_off_stop  = lights_off_stops[i]

                pat.data.dome_flat_field(flat_date, band, sftp=sftp, upload=True, delete_raw=True, lights_on_start=lights_on_start, lights_on_stop=lights_on_stop, lights_off_start=lights_off_start, lights_off_stop=lights_off_stop)

    print('')

    #Variable pixel masks. 
    if variable_pixels:
        for i in range(len(exptimes)):
            dark_date = dark_dates[i]
            exptime = float(exptimes[i])
            vpm_path = calibration_path/('Variable Pixel Masks/vpm_'+str(exptime)+'_s_'+dark_date+'.fits')
            if vpm_path.exists():
                print('{} already exists in {}, skipping.'.format(vpm_path.name, calibration_path))
            else:
                print('Making variable pixel mask for {}-s exposure time on {}.'.format(exptime, dark_date))
                pat.data.variable_pixels(dark_date, exptime, upload=True, sftp=sftp)

    print('')

    #Hot pixel masks. 
    if hot_pixels:
        for i in range(len(exptimes)):
            dark_date = dark_dates[i]
            exptime = float(exptimes[i])
            hpm_path = calibration_path/('Hot Pixel Masks/hpm_'+str(exptime)+'_s_'+dark_date+'.fits')
            if hpm_path.exists():
                print('{} already exists in {}, skipping.'.format(hpm_path.name, calibration_path))
            else:
                print('Making hot pixel mask for {}-s exposure time on {}.'.format(exptime, dark_date))
                pat.data.hot_pixels(dark_date, exptime, box_l=5, upload=True, sftp=sftp)

    print('')
    
    #Dead pixel masks
    if dead_pixels:
        for i in range(len(bands)):
            flat_date = flat_dates[i]
            band = bands[i]
            dpm_path = calibration_path/('Dead Pixel Masks/dpm_'+band+'_'+flat_date+'.fits')
            if dpm_path.exists():
                print('{} already exists in {}, skipping.'.format(dpm_path.name, calibration_path))
            else:
                print('Making dead pixel mask for {}-band exposure time on {}.'.format(band, dark_date))
                pat.data.dead_pixels(flat_date, band, upload=True, sftp=sftp)
    
    print('')
    
    #Bad pixel masks. 
    if bpm:
        for i in range(len(exptimes)):
            dark_date = dark_dates[i]
            exptime = float(exptimes[i])
            for j in range(len(bands)):
                flat_date = flat_dates[j]
                band = bands[j]
                bpm_path = calibration_path/('Bad Pixel Masks/bpm_'+band+'_'+str(exptime)+'_s_'+flat_date+'.fits')
                if bpm_path.exists():
                    print('{} already exists in {}, skipping.'.format(bpm_path.name, calibration_path))
                else:
                    print('Making dead pixel mask for {}-band, {}-s exposure time on {}.'.format(band, exptime, flat_date))
                    pat.data.bpm_maker(flat_date, dark_date, exptime, band, upload=True, sftp=sftp)
    return 