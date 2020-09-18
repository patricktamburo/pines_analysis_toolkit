import pines_analysis_toolkit as pat

date = '20200907'
band = 'J'
exptime = 60.0

sftp = pat.utils.pines_login()
pat.data.bpm_maker(date, exptime, band, upload=True, sftp=sftp)
sftp.close()