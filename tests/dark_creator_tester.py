import pines_analysis_toolkit as pat
import pdb

date = '20200907'
exptime = 60.
sftp = pat.utils.pines_login()
pat.data.dark(sftp, date, exptime, upload=True, delete_raw=True)
sftp.close()