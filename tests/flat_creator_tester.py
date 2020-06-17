import pines_analysis_toolkit as pat

date = '20200512'
band = 'H'
pat.data.dome_flat_field(date, band, lights_on_start=1, lights_on_stop=100, lights_off_start=101, lights_off_stop=200, upload=True, delete_raw=True)