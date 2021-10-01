
def get_calibrations(sftp, pines_path):
    """Downloads all calibration files on the PINES server. 

    :param sftp: sftp connection to the PINES server.
    :type sftp: pysftp connection
    :param pines_path: path to top-level PINES_analysis_toolkit directory containing the 'Calibrations/' folder
    :type pines_path: pathlib.PosixPath
    """

    print('Downloading calibration files from pines.bu.edu!')
    print('')
    
    sftp.chdir('/')
    sftp.chdir('data/calibrations/')
    calibration_directories = ['Bad Pixel Masks', 'Darks', 'Dead Pixel Masks', 'Flats', 'Hot Pixel Masks', 'Kokopelli Mask', 'Variable Pixel Masks']
    for i in range(len(calibration_directories)):
        sftp.chdir(calibration_directories[i])
        if calibration_directories[i] != 'Flats':
            files = sftp.listdir()
            #Set the directory the calibrations will be downloaded to.
            if calibration_directories[i] == 'Darks':
                local_path = pines_path/('Calibrations/'+calibration_directories[i]+'/Master Darks/')
            else:
                local_path = pines_path/('Calibrations/'+calibration_directories[i]+'/')

            for j in range(len(files)):
                if not (local_path/files[j]).exists():
                    sftp.get(files[j], local_path/files[j])
                    print('Downloading to {}'.format(local_path/files[j]))
        else:
            sftp.chdir('Domeflats')
            bands = sftp.listdir()
            for j in range(len(bands)):
                sftp.chdir(bands[j])
                files = sftp.listdir()
                local_path = pines_path/('Calibrations/Flats/Domeflats/'+bands[j]+'/Master Flats/')
                for k in range(len(files)):
                    if not (local_path/files[k]).exists():
                        sftp.get(files[k], local_path/files[k])
                        print('Downloading to {}'.format(local_path/files[k]))
                sftp.chdir('..')
        sftp.chdir('/')
        sftp.chdir('data/calibrations/')

