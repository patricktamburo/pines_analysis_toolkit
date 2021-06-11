import pines_analysis_toolkit as pat 

sftp = pat.utils.pines_login() 

targets = ['2MASS J15394189-0520428'
'2MASS J15551573-0956055',
'2MASS J15552614+0017204',
'2MASS J16051015+0421521',
'2MASS J16192830+0050118']

#for target in targets:
#  pat.data.get_reduced_science_files(sftp, target)

#Declare list of all UT OBSERVING DATES for which you want to update logs. 
dates = ['20210527', '20210528'] 

#Update logs with more accurate shifts.
for date in dates:
   pat.observing.log_updater(date, sftp, upload=False)

pat.data.straggler_upload(sftp, targets)
sftp.close()