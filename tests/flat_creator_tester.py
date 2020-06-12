import pines_analysis_toolkit as pat

run = '20200307'
band = 'H'
pat.data.dome_flat_field(run, band, delete_raw=True, upload=True)