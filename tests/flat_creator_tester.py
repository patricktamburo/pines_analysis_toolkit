import pines_analysis_toolkit as pat

date = '20200907'
band = 'J'
sftp = pat.utils.pines_login()
pat.data.dome_flat_field(sftp, date, band, upload=True, delete_raw=True)
sftp.close()