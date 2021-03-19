import pines_analysis_toolkit as pat 

targets = ['2MASS J11394192-0310039', '2MASS J11491231-0153003',
           '2MASS J11555389+0559577', '2MASS J11582077+0435014']

sftp = pat.utils.pines_login()

for target in targets:
    pat.data.get_raw_science_files(sftp, target)
    pat.data.reduce(target, delete_raw=True, delete_reduced=False, upload=True, sftp=sftp)

sftp.close()

