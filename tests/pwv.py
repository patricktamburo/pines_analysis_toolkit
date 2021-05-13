from pines_analysis_toolkit.pwv.fyodor import download_nc, rename_nc
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
import pdb 

ut_dates = ['20201203', '20201204', '20201205', '20201206']
url_codes = ['8042323663', '8042323665', '8042323666', '8042323878']

pines_path = pines_dir_check()
download_directory = pines_path/('Calibrations/PWV/')

if len(ut_dates) != len(url_codes):
    raise RuntimeError('There must be one url code for every UT date you want to download.')

for i in range(len(ut_dates)):
    date = ut_dates[i]
    url = 'https://download.avl.class.noaa.gov/download/'+url_codes[i]+'/001'
    download_nc(url, str(download_directory), date, n_threads=1)
    rename_nc(str(download_directory/date))

