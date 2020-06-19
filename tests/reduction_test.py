import pines_analysis_toolkit as pat
import pdb

target = '2MASS J13334540-0215599'
sftp = pat.utils.pines_login()
pat.data.get_raw_science_files(sftp, target)
pat.data.reduce(target, upload=True, sftp=sftp, delete_raw=True, delete_reduced=False)
sftp.close()
