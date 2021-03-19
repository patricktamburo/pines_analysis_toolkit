import pines_analysis_toolkit as pat 

#Declare list of all UT OBSERVING DATES for which you want to update logs. 
dates = ['20200303', '20200304', '20200305', '20200306', '20200307']

sftp = pat.utils.pines_login() 

#Update logs with more accurate shifts.
for date in dates:
    pat.observing.log_updater(date, sftp, upload=True)

sftp.close()