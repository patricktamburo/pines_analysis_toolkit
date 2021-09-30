import pines_analysis_toolkit as pat

date = '20200907'
band = 'J'
sftp = pat.utils.pines_login()
#pat.data.dome_flat_field(date, band, upload=True, sftp=sftp)
pat.data.dead_pixels(date, band, upload=True, sftp=sftp)
sftp.close()
