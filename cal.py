import pines_analysis_toolkit as pat 
import pdb 

sftp = pat.utils.pines_login()

#Set up calibration info. In this example, I'll make calibrations for the 202005 run, which had targets observed in H band with 10- and 15-s exposure times. Get calibration info from PATATO: https://docs.google.com/spreadsheets/d/1ACMGHrXuVh6yyzZWsR-dy3bq72bC0ElfxwUnfeA8HDM/edit#gid=0. 

#Dark info. Each exposure time must have a corresponding dark date!
exptimes =   [        15]
dark_dates = ['20200307']

#Flat info. Each band must have a corresponding flat date! 
bands      = [       'H']
flat_dates = ['20200307']
lights_on_starts  = [386]
lights_on_stops   = [585]
lights_off_starts = [586]
lights_off_stops  = [785] 

#Run the calibrations. NOTE that if a calibration file already exists on your local computer, it will NOT create a new copy of it. If you want to make a new copy, you will need to remove the existing calibration file from the appropriate ~/Documents/PINES_analysis_toolkit/Calibrations/ directory. Also NOTE that any created master calibrations will automatically be uploaded to the PINES server! 
pat.data.make_calibrations(sftp, exptimes, bands, dark_dates, flat_dates, dark_starts=[], dark_stops=[], lights_on_starts=lights_on_starts, lights_on_stops=lights_on_stops, lights_off_starts=lights_off_starts, lights_off_stops=lights_off_stops)

sftp.close()