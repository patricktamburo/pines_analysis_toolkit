import pines_analysis_toolkit as pat

run = '20200531'
band = 'J'
pat.data.dome_flat_field(run, band, delete_raw=True)