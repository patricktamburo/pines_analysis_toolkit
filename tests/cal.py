import pines_analysis_toolkit as pat 
import pdb 

sftp = pat.utils.pines_login()

#Set up calibration info. In this example, I'll make calibrations for the 202005 run, which had targets observed in H band with 10- and 15-s exposure times. Get calibration info from PATATO: https://docs.google.com/spreadsheets/d/1ACMGHrXuVh6yyzZWsR-dy3bq72bC0ElfxwUnfeA8HDM/edit#gid=0. 

#Dark info. Each exposure time must have a corresponding dark date!
exptimes =   [        10,         15]
dark_dates = ['20200512', '20200512']

#Flat info. Each band must have a corresponding flat date! NOTE that despite there being multiple exposure times, all of the targets were observed in the same band (H), so we only need to make one flat. 
bands      = ['H']
flat_dates = ['20200512']

#Run the calibrations. NOTE that if a calibration file already exists on your local computer, it will NOT create a new copy of it. If you want to make a new copy, you will need to remove the existing calibration file from the appropriate ~/Documents/PINES_analysis_toolkit/Calibrations/ directory. Also NOTE that any created master calibrations will automatically be uploaded to the PINES server! 
pat.data.make_calibrations(sftp, exptimes, bands, dark_dates, flat_dates)

sftp.close()