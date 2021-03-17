import pines_analysis_toolkit as pat 
import pdb 

def make_calibrations(sftp, exptimes, bands, dark_dates, flat_dates, dark=True, flat=True, variable_pixels=True, 
                      hot_pixels=True, dead_pixels=True, bpm=True):
    '''Authors: 
        Patrick Tamburo, Boston University, March 2021
    Purpose: 
        Makes all calibration files for a target. 
    Inputs:
        sftp (pysftp object): connection to the PINES server. 
        exptimes (list of floats): list of all exposure times for this target (e.g., [15., 30.] if observed with both 15 
            and 30 second exposure times).
        bands (list of str): list of all bands for this target (e.g., ['J', 'H'] if observed in both J and H bands).
        dark_dates (list of str): list of UT dates on which darks were taken (order matches order of exptimes!).
        flat_dates (list of str): list of UT dates on which flats were taken (order matches order of bands!)
        dark (bool): Whether or not to make the master darks. 
        flat (bool): Whether or not to make the master flats. 
        variable_pixels (bool): Whether or not to make the variable pixel masks. 
        hot_pixels (bool): Whether or not to make the hot pixel masks. 
        dead_pixels (bool): Whether or not to make the dead pixel masks. 
        bpm (bool): Whether or not to make the bad pixel masks.  
    Outputs:
        Writes calibration images to appropriate ~/Documents/PINES_analysis_toolkit/Calibrations/ directory. 
    TODO:
        Have to allow for dark_start/stop, lights_on_start/stop, lights_off_start/stop arguments in case you need to specify numbers of dark/flat files. 
    '''

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
                pat.data.dark(dark_date, exptime, upload=True, delete_raw=True, sftp=sftp, dark_start=0, dark_stop=0)

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
                pat.data.dome_flat_field(flat_date, band, sftp=sftp, upload=True, delete_raw=True, lights_on_start=0, lights_on_stop=0, lights_off_start=0, lights_off_stop=0)

    print('')

    #Variable pixel masks. 
    if variable_pixels:
        for i in range(len(exptimes)):
            dark_date = dark_dates[i]
            exptime = float(exptimes[i])
            vpm_path = calibration_path/('Variable Pixel Masks/vpm'+str(exptime)+'_s_'+dark_date+'.fits')
            if vpm_path.exists():
                print('{} already exists in {}, skipping.'.format(vpm_path.name, calibration_path))
            else:
                print('Making variable pixel mask for {}-s exposure time on {}.'.format(exptime, dark_date))
                pat.data.variable_pixels(dark_date, exptime, upload=True, sftp=sftp)

    #Hot pixel masks. 
    if hot_pixels:
        for i in range(len(exptimes)):
            dark_date = dark_dates[i]
            exptime = float(exptimes[i])
            hpm_path = calibration_path/('Hot Pixel Masks/hpm'+str(exptime)+'_s_'+dark_date+'.fits')
            if hpm_path.exists():
                print('{} already exists in {}, skipping.'.format(hpm_path.name, calibration_path))
            else:
                print('Making hot pixel mask for {}-s exposure time on {}.'.format(exptime, dark_date))
                pat.data.hot_pixels(dark_date, exptime, box_l=5, upload=True, sftp=sftp)

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
    
    #Bad pixel masks. 
    if bpm:
        for i in range(len(exptimes)):
            dark_date = dark_dates[i]
            exptime = exptimes[i]
            for j in range(len(bands)):
                flat_date = flat_dates[j]
                band = bands[j]
                bpm_path = calibration_path/('Bad Pixel Masks/bpm_'+band+'_'+str(exptime)+'_s_'+flat_date+'.fits')
                if bpm.path.exists():
                    print('{} already exists in {}, skipping.'.format(bpm_path.name, calibration_path))
                else:
                    print('Making dead pixel mask for {}-band, {}-s exposure time on {}.'.format(band, exptime, flat_date))
                    pat.data.bpm_maker(flat_date, dark_date, exptime, band, upload=True, sftp=sftp)
    return 