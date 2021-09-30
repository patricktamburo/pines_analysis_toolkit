from pines_analysis_toolkit.pwv.fyodor import download_nc, rename_nc
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
import pdb 

ut_dates = ['20200628', '20200630', '20200701', '20200702']
url_codes = ['8066550611','8066550612','8066550613','8066550614']

pines_path = pines_dir_check()
download_directory = pines_path/('Calibrations/PWV/')

if len(ut_dates) != len(url_codes):
    raise RuntimeError('There must be one url code for every UT date you want to download.')

for i in range(len(ut_dates)):
    date = ut_dates[i]
    url = 'https://download.avl.class.noaa.gov/download/'+url_codes[i]+'/001'
    download_nc(url, str(download_directory), date, n_threads=1)
    rename_nc(str(download_directory/date))

