import pines_analysis_toolkit as pat 
import pandas as pd 
from datetime import datetime
import julian 
from scipy import interpolate
import numpy as np 
import pdb 
import matplotlib.pyplot as plt 
plt.ion() 

def los_pwv_data_generator(target, ut_dates, centroided_sources, interp_type='linear'):
    pines_path = pat.utils.pines_dir_check()
    short_name = pat.utils.short_name_creator(target)
    pines_sample_path = pines_path/('Misc/PINES sample.xlsx')
    sample_df = pd.read_excel(pines_sample_path)
    target_row = np.where(np.array(sample_df['2MASS Name']) == target.split('J')[1])[0][0]
    target_ra = sample_df['RA (deg)'][target_row]
    target_dec = sample_df['Dec (deg)'][target_row]
    for date in ut_dates:
        fyodor_data_path = pines_path/('Calibrations/PWV/'+date)
        pat.pwv.fyodor.pwv(target, fyodor_data_path, P_min=785, P_max=300, RA=target_ra, Dec=target_dec, csv=True)
    
    #Interpolate PWV data onto grid of PINES times
    time_strs = []
    pwv   = []
    for date in ut_dates:
        fyodor_data_path = pines_path/('Objects/'+short_name+'/pwv/')
        csv_path = fyodor_data_path/('PWV_los_'+date+'.csv')
        df = pd.read_csv(csv_path)
        time_strs.extend(df['Time'])
        pwv.extend(df['PWV'])
    

    time_strs = np.array(time_strs)
    fyodor_dates = np.array([datetime.strptime(i, '%Y-%m-%d %H:%M:%S') for i in time_strs]) 
    fyodor_times = np.array([julian.to_jd(i) for i in fyodor_dates]) #Convert to JD UTC
    pwv = np.array(pwv)
    fyodor_times = fyodor_times[~np.isnan(pwv)]
    pwv = pwv[~np.isnan(pwv)]

    centroided_sources.columns = centroided_sources.columns.str.strip()
    pines_times = centroided_sources['Time (JD UTC)']

    #Interpolate Fyodor full disk pwv data onto the PINES data grid.
    if interp_type == 'linear':
        f1 = interpolate.interp1d(fyodor_times, pwv, kind='nearest')
        pwv_interp = f1(pines_times)

    elif interp_type == 'cubicspline':
        tck = interpolate.splrep(fyodor_times, pwv)
        pwv_interp = interpolate.splev(pines_times, tck, der=0)

    pwv_output_dict = {'Time (JD UTC)':pines_times, 'PWV':pwv_interp}
    pwv_output_df = pd.DataFrame(pwv_output_dict)
    pwv_output_path = pines_path/('Objects/'+short_name+'/pwv/'+short_name+'_fyodor_pwv.csv')
    pwv_output_df.to_csv(pwv_output_path, index=0)
    