import pines_analysis_toolkit as pat
import pdb

date = '20200907'
exptime = 60.
sftp = pat.utils.pines_login()
#pat.data.dark(date, exptime, upload=True, delete_raw=False, sftp=sftp)
#pat.data.variable_pixels(date, exptime, upload=True, sftp=sftp)
pat.data.hot_pixels(date, exptime, upload=True, sftp=sftp)
sftp.close()
