import pines_analysis_toolkit as pat
import pdb

target = '2MASS J23255604-0259508'
#sftp = pat.utils.pines_login()
#pat.data.get_raw_science_files(sftp, target)
pat.data.reduce(target, upload=False, delete_raw=False, delete_reduced=False)
#sftp.close()
