import pines_analysis_toolkit as pat

date = '20200531'
exptime = 60.
pat.data.dark(date, exptime, upload=True, delete_raw=True)