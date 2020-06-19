import pines_analysis_toolkit as pat
import pdb

date = '20200512'
exptime = 10.
sftp = pat.utils.pines_login()
pat.data.dark(sftp, date, exptime, dark_start=221, dark_stop=240, upload=True, delete_raw=True)
sftp.close()