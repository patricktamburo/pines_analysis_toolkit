import pines_analysis_toolkit as pat 

sftp = pat.utils.pines_login() 

# targets = ['2MASS J08312221+1538511',
# '2MASS J08350622+1953050',
# '2MASS J08354537+2224310',
# '2MASS J08433323+1024470',
# '2MASS J08560211+1240150']

#for target in targets:
#  pat.data.get_reduced_science_files(sftp, target)

#Declare list of all UT OBSERVING DATES for which you want to update logs. 
dates = ['20210527', '20210528'] 

#Update logs with more accurate shifts.
for date in dates:
   pat.observing.log_updater(date, sftp, upload=False)

#pat.data.straggler_upload(sftp, targets)
sftp.close()