
def get_master_log(sftp, pines_path):
    """Downloads master_log.txt from the PINES server.

    :param sftp: sftp connection to the PINES server
    :type sftp: pysftp connection
    :param pines_path: top-level PINES path
    :type pines_path: pathlib.PosixPath
    """

    print('Downloading master_log.txt to {}/Logs/'.format(pines_path))
    sftp.chdir('/code/PINES_server_scripts/')
    sftp.get('master_log.txt',(pines_path/'Logs/master_log.txt'))
    return