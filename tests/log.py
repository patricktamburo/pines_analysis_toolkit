import pines_analysis_toolkit as pat 

sftp = pat.utils.pines_login() 

targets = ['2MASS J17281134+0839590'] #13

for target in targets:
  pat.data.get_reduced_science_files(sftp, target)

#Declare list of all UT OBSERVING DATES for which you want to update logs. 
dates = ['20200716'] 

#Update logs with more accurate shifts.
for date in dates:
   pat.observing.log_updater(date, sftp, upload=False)

pat.data.straggler_upload(sftp, targets)
sftp.close()