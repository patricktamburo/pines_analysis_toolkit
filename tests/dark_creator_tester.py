import pines_analysis_toolkit as pat

date = '20200512'
exptime = 10.
pat.data.dark(date, exptime, dark_start=221, dark_stop=240, upload=True, delete_raw=True)