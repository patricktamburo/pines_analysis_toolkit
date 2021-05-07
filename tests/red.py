import pines_analysis_toolkit as pat 

targets = ['2MASS J23062928-0502285']

sftp = pat.utils.pines_login()

for target in targets:
    pat.data.get_raw_science_files(sftp, target)
    pat.data.reduce(target, delete_raw=True, delete_reduced=False, upload=True, sftp=sftp)

sftp.close()

 