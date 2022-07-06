from pines_analysis_toolkit.utils import pines_dir_check, short_name_creator, jd_utc_to_bjd_tdb, pines_log_reader

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import EarthLocation, AltAz
from astropy.time import Time
from astropy.utils.data import clear_download_cache
from astropy.constants import R_earth
from astropy.io import fits
from astropy.convolution import convolve, Gaussian1DKernel

from netCDF4 import Dataset
import os
import numpy as np
import time
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import requests
from bs4 import BeautifulSoup
from pines_analysis_toolkit.utils import get_source_names
import wget
from multiprocessing.pool import ThreadPool
from datetime import datetime
import julian
from glob import glob
import re
import math
from pathlib import Path

from scipy import interpolate
import scipy.signal
from scipy.ndimage.filters import maximum_filter, uniform_filter
from scipy.stats import norm 
   
#__all__ = ['rename_nc', 'download_nc', 'pwv']

clear_download_cache()


def rename_nc(directory):
    """
    Rename GOES-R .nc files downloaded via subscription. It removes the order number at the start of the filename. 
    This must be done so that the pwv function can retrieve measurements in chronological order.

    Parameters
    ----------
    directory : str
        Directory path containing the files.
    """
    full_direc = os.listdir(directory)
    for i in range(len(full_direc)):
        try:
            tmp = int(full_direc[i].split('.')[0])
            dst = full_direc[i].split('.')[1] + '.' + \
                full_direc[i].split('.')[2]
            src = full_direc[i]
            os.rename(src, dst)
        except Exception:
            pass


def create_folder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
        else:
            os.chdir(directory)
    except OSError:
        print("Error: Creating directory." + directory)


def download_file_wget(url):
    local_filename = url.split('/')[-1]
    wget.download(url)
    return url


def download_nc(url, directory, date, n_threads=8):
    """
    Download .nc files from the NOAA webpage after manually requesting the dataset and store them in a new folder 

    Parameters
    ----------
    url : str
        Website provided by NOAA with the requested files
    directory : str
        Path where a folder with the files should be created
    date: str, Format yyyymmdd
        Day when the measurement took place
    n_threads: int, number of threads to use to download files
    """

    any_day = str(date)
    date = time.strptime(any_day, '%Y%m%d')

    create_folder(str(directory) + "/{}".format(any_day))
    os.chdir(directory + '/{}'.format(any_day))

    url = str(url)
    r = requests.get(url, allow_redirects=False)
    soup = BeautifulSoup(r.content, 'html.parser')
    TableContent = soup.select('table tr a')

    Content = []
    for i in range(0, len(TableContent)):
        Contentt = TableContent[i]['href'].split('/')[-1]
        Content.append(Contentt)

#    if len(str(date[7])) == 1:
#        subs = '{}'.format(date[0]) + '00' + '{}'.format(date[7])
#    elif len(str(date[7])) == 2:
#        subs = '{}'.format(date[0]) + '0' + '{}'.format(date[7])
#    elif len(str(date[7])) == 3:
#        subs = '{}'.format(date[0]) + '{}'.format(date[7])

    # Modified this so that we know how many files to download before we start.
    FileName = [i for i in Content]
    FileName = [path for path in FileName if not os.path.exists('./' + path)]

    urls = []
    for i in range(0, len(FileName)):
        urlstemp = url + "/" + FileName[i]
        urls.append(urlstemp)
    
    print('There are %s files to download' % len(urls))

    for ii in range(0, len(urls), n_threads):
        slice_item = slice(ii, ii + n_threads, 1)

        files = ThreadPool(len(urls[slice_item])).imap_unordered(
            download_file_wget, urls[slice_item])
        # Need to do this print so that it waits before starting the next batch.
        for f in files:
            print('%s complete' % f.split('/')[-1])

    print("All files downloaded!")


def fyodor_pwv(short_name, directory, P_min, P_max, line_of_sight='target', RA=None, Dec=None, plot=False, csv=False, force_output_path=''):
    """
    Compute the precipitable water vapor (PWV) at the PTO at zenith or in direction of
    ``target``.

    Parameters
    ----------
    directory : str
        Working directory with GOES-R files in it
    P_min : float
        Lower pressure level (lower altitude). Range between 1100 and 0.05 (hPa)
    P_max : float
        Upper pressure level (higher altitude). Range between 1100 and 0.05 (hPa)
    line_of_sight :  str, {"target", "zenith"}
        Either compute line of sight to the target or to the zenith.
    RA : float
        Right ascension of target (in degrees)
    Dec : float
        Declination of target (in degrees)
    plot : bool
        Generate a plot of the PWV at each time.
    csv : bool
        Generate a csv file of the PWV at each time. 

    Returns
    -------
    dates : list
    PWV : `~numpy.ndarray`
    LVT : `~numpy.ndarray`
    LVM : `~numpy.ndarray`
    """

    # Coordinates input in degrees
    latdeg = 35.097289793027436
    londeg = -111.53686622502691
    site = 'Perkins Telescope Observatory'
    
    #Fix the names of files so they appear in order. 
    rename_nc(directory)

    # Access working directory

    # nc_filesT = glob.glob('*OR_ABI-L2-LVTPF*') #MODIFIED TO WORK WITH CONUS DATA
    nc_filesT = glob(str((directory/'*OR_ABI-L2-LVTP*.nc')))

    nc_filesT = sorted(nc_filesT)
    #nc_filesM = glob.glob('*OR_ABI-L2-LVMPF*')
    nc_filesM = glob(str((directory/'*OR_ABI-L2-LVMP*.nc')))
    nc_filesM = sorted(nc_filesM)
    
    # Open the first file to retrieve earth parameters
    Proj_info = Dataset(nc_filesT[0], 'r')
    proj_info = Proj_info.variables['goes_imager_projection']
    lon_origin = proj_info.longitude_of_projection_origin
    H = proj_info.perspective_point_height+proj_info.semi_major_axis
    r_eq = proj_info.semi_major_axis
    r_pol = proj_info.semi_minor_axis

    # Retrieve pressure data
    P = Proj_info.variables['pressure'][:]
    Proj_info.close()
    Proj_info = None

    # Retrieve time data
    g16_data_file = []
    t = []
    epoch = []
    date = []
    day = []
    hour = []

    for i in range(0, len(nc_filesT)):
        g16_data = nc_filesT[i]
        g16_data_file.append(g16_data)
        g16 = Dataset(g16_data_file[i], 'r')
        ttemp = g16.variables['t'][:]
        t.append(ttemp)
        epochtemp = 946728000 + int(t[i])
        epoch.append(epochtemp)
        datetemp = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(epoch[i]))
        date.append(datetemp)
        daytemp = time.strftime("%d-%m-%Y", time.gmtime(epoch[i]))
        day.append(daytemp)
        hourtemp = time.strftime("%H:%M:%S", time.gmtime(epoch[i]))
        hour.append(hourtemp)

    # Get the user's pines_analysis_toolkit path
    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()

    out_date = day[1][-4:]+day[1][3:5]+day[1][0:2]
    # Make the object's pwv directory if it doesn't already exist.
    if not os.path.isdir(pines_path/('Objects/'+short_name+'/pwv/')):
        os.mkdir(pines_path/('Objects/'+short_name+'/pwv/'))
    output_path = pines_path / \
        ('Objects/'+short_name+'/pwv/PWV_los_{}.csv'.format(out_date))
    if os.path.exists(output_path):
        # print('PWV output already exists for {}, returning.'.format(out_date))
        print('')
        #return

    # Use astropy.time to keep format for target coordinates:
    times = Time(date, format='iso', scale='utc')

    # Barometric formula
    p0 = P[0]  # hPa
    R_D = 287  # Jkg-1K-1
    g = 9.81  # m/s2
    T_s = 288  # K
    L = -0.0065  # K/m

    h = (T_s/L)*((P/p0)**(-L*R_D/g)-1)*u.m

    e = np.sqrt((r_eq**2-r_pol**2)/(r_eq**2))  # Eccentricity

    latdeg = float(latdeg)
    londeg = float(londeg)

    # Pressure level boundaries
    P_minb = np.abs(P-P_min).argmin()
    P_maxb = np.abs(P-P_max).argmin()

    loc = EarthLocation(lat=latdeg*u.degree, lon=londeg*u.degree)

    # Convert from radian to degrees:
    raddeg = 180/np.pi

    if line_of_sight == 'zenith':
        latt = latdeg
        lont = londeg
    elif line_of_sight == 'target':
        RA = float(RA)
        Dec = float(Dec)
        Sky = SkyCoord(ra=RA*u.degree, dec=Dec*u.degree)
        Aa = Sky.transform_to(AltAz(obstime=times, location=loc))
        latt = Dec
        lont = RA

    # Computes PWV along line of sight
    if line_of_sight == 'target':
        # ORIGINAL USED VALUES WERE WHERE ALT WAS >30. CHANGED FOR PERKINS.
        INDEX = np.ravel(np.where(Aa.alt.degree < 0))
        INDEXP = np.ravel(np.where(Aa.alt.degree > 0))

        # Keep time values corresponding to Alt above 30 degrees
        EPOCH = epoch
        for index in sorted(INDEX, reverse=True):
            del EPOCH[index]
        HOUR = []
        for i in range(0, len(epoch)):
            HOURt = time.strftime("%H:%M:%S", time.gmtime(EPOCH[i]))
            HOUR.append(HOURt)
        DATE = []
        for i in range(0, len(epoch)):
            DATEt = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(EPOCH[i]))
            DATE.append(DATEt)

        Alt = Aa.alt.rad
        Az = Aa.az.rad

        # Compute distance from location to projection point
        delta_x = []
        d_lat = []
        d_lon = []

        for i in range(0, len(Alt)):
            delta_xi = []
            for j in h:
                delta_xt = j/np.tan(Alt[i])
                delta_xt = delta_xt*u.m**-1
                delta_xi.append(delta_xt)
            delta_x.append(delta_xi)
        delta_x = np.array(delta_x)

        for i in range(0, len(Az)):
            d_latt = delta_x[i, ]*np.cos(Az[i])
            d_lont = delta_x[i, ]*np.sin(Az[i])
            d_lat.append(d_latt)
            d_lon.append(d_lont)
        d_lat = np.array(d_lat)
        d_lon = np.array(d_lon)

        # Compute latitude and longitude of projection points
        lat_proj = []
        lon_proj = []
        for i in range(0, len(Alt)):
            obs_latt = loc.lat.degree + raddeg * \
                (np.arctan(d_lat[i, ]/R_earth*u.m**1)*u.rad**-1)
            obs_lont = loc.lon.degree + raddeg * \
                (np.arctan(d_lon[i, ]/R_earth*u.m**1)*u.rad**-1)
            lat_proj.append(obs_latt)
            lon_proj.append(obs_lont)
        lat_proj = np.array(lat_proj)
        lon_proj = np.array(lon_proj)
        rad = (np.pi)/180

        lat_proj_rad = rad*lat_proj
        lon_proj_rad = rad*lon_proj
        lambda_0 = rad*lon_origin

        # T ransform into scan angles
        lat_origin = np.arctan(((r_pol**2)/(r_eq**2))*np.tan(lat_proj_rad))

        r_c = r_pol/(np.sqrt(1-(e**2)*(np.cos(lat_origin))**2))

        s_x = H - r_c*np.cos(lat_origin)*np.cos(lon_proj_rad-lambda_0)
        s_y = -r_c*np.cos(lat_origin)*np.sin(lon_proj_rad-lambda_0)
        s_z = r_c*np.sin(lat_origin)

        s = np.sqrt(s_x**2+s_y**2+s_z**2)

        x = np.arcsin(-s_y/s)
        y = np.arctan(s_z/s_x)

        g16_data_fileT = []
        xscanT = []
        yscanT = []
        XT = []
        YT = []
        LVT = []

        # Retrieve Temperature data
        for i in INDEXP:
            g16_data = nc_filesT[i]
            g16_data_fileT.append(g16_data)
            g16 = Dataset(nc_filesT[i], 'r')
            xtemp = g16.variables['x'][:]
            xscanT.append(xtemp)
            ytemp = g16.variables['y'][:]
            yscanT.append(ytemp)

            LVTi = []
            Xi = []
            Yi = []
            for j in range(P_minb, P_maxb+1):
                Xtemp = np.abs(xtemp-x[i, j]).argmin()
                Xi.append(Xtemp)
                Ytemp = np.abs(ytemp-y[i, j]).argmin()
                Yi.append(Ytemp)
                LVTtemp = g16.variables['LVT'][Ytemp, Xtemp, j]
                LVTi.append(LVTtemp)
            LVT.append(LVTi)
            XT.append(Xi)
            YT.append(Yi)
        LVT = np.array(LVT)

        # Retrieve Relative humidity data
        g16_data_fileM = []
        xscanM = []
        yscanM = []
        XM = []
        YM = []
        LVM = []

        for i in INDEXP:
            g16_dataM = nc_filesM[i]
            g16_data_fileM.append(g16_dataM)
            g16M = Dataset(nc_filesM[i], 'r')
            xtempM = g16M.variables['x'][:]
            xscanM.append(xtempM)
            ytempM = g16M.variables['y'][:]
            yscanM.append(ytempM)

            LVMi = []
            Xi = []
            Yi = []
            for j in range(P_minb, P_maxb+1):
                XtempM = np.abs(xtempM-x[i, j]).argmin()
                Xi.append(XtempM)
                YtempM = np.abs(ytempM-y[i, j]).argmin()
                Yi.append(YtempM)
                LVMtemp = g16M.variables['LVM'][YtempM, XtempM, j]
                LVMi.append(LVMtemp)
            LVM.append(LVMi)
            XM.append(Xi)
            YM.append(Yi)
        LVM = np.array(LVM)

        P = 100*P
        LVT = LVT-273.15

        Pi = P[P_minb:P_maxb+1]

        # Constants needed and integrand
        rho_w = 1000  # kg/m**3
        g = 9.81  # m/s**2
        C = (-1)/(rho_w*g)
        # Partial water vapour pressure in Pa
        ev = 100*6.112*LVM*np.exp((17.67*LVT)/(LVT+243.5))
        q = (0.622*ev)/(Pi-0.378*ev)  # Specific humdity
        # Complete integrand multiplied by 1000 to get the PWV in mm.
        f = 1000*C*q

        # Numerical integration
        PWV = []
        for j in range(0, len(LVT)):
            integral = 0
            for i in range(1, len(Pi)):
                integral = integral + (Pi[i]-Pi[i-1])*((f[j, i]+f[j, i-1])/2)
            PWV.append(integral)

        PWV = np.asarray(PWV)

        if plot:
            # Plot and save data
            fig = plt.figure(figsize=(15, 10))
            ax = fig.add_subplot(111)
            ax.plot(HOUR, PWV, 'bo', ms=4)
            plt.title('PWV along line of sight, {} on {}'.format(
                site, day[1]), fontsize=26)
            plt.xticks(rotation='vertical', fontsize=24)
            plt.yticks(fontsize=24)
            ax.set_xlabel("Date", color="C0", fontsize=24)
            ax.set_ylabel("PWV (mm)", color="C0", fontsize=24)
            RA_patch = mpatches.Patch(
                color='white', label='RA: {} degrees'.format(RA))
            Dec_patch = mpatches.Patch(
                color='white', label='Dec: {} degrees'.format(Dec))
            every_nth = 6
            for n, label in enumerate(ax.xaxis.get_ticklabels()):
                if n % every_nth != 0:
                    label.set_visible(False)
            for n, label in enumerate(ax.xaxis.get_ticklines()):
                if n % every_nth != 0:
                    label.set_visible(False)
            plt.tight_layout()
            plt.legend(handles=[RA_patch, Dec_patch],
                       loc='lower right', fontsize=22)
            plt.show()
            fig.savefig('PWV_line_of_sight_{}_{}.png'.format(site, day[1]))

        if csv:
            np.savetxt(output_path, np.column_stack((date, PWV)),
                       delimiter=',', fmt='%s', header='Time,PWV', comments='')

    # Computes PWV at zenith
    elif line_of_sight == 'zenith':

        # Transform latitude and longitude into scan angles
        rad = np.pi/180
        lambda_0 = rad*lon_origin
        obs_lat_rad = rad*latt
        obs_lon_rad = rad*lont

        lat_origin = np.arctan(((r_pol**2)/(r_eq**2))*np.tan(obs_lat_rad))

        r_c = r_pol/(np.sqrt(1-(e**2)*(np.cos(lat_origin))**2))

        s_x = H - r_c*np.cos(lat_origin)*np.cos(obs_lon_rad-lambda_0)
        s_y = -r_c*np.cos(lat_origin)*np.sin(obs_lon_rad-lambda_0)
        s_z = r_c*np.sin(lat_origin)

        s = np.sqrt(s_x**2+s_y**2+s_z**2)

        x = np.arcsin(-s_y/s)
        y = np.arctan(s_z/s_x)

        xscanT = []
        yscanT = []

        # Retrieve Temperature data
        LVT = []
        for i in range(0, len(nc_filesT)):
            g16_data = nc_filesT[i]
            g16_data_file.append(g16_data)
            g16 = Dataset(g16_data_file[i], 'r')
            xtemp = g16.variables['x'][:]
            xscanT.append(xtemp)
            ytemp = g16.variables['y'][:]
            yscanT.append(ytemp)

            XT = []
            YT = []
            LVTi = []
            for j in range(0, len(P)):
                Xtemp = np.abs(xtemp-x).argmin()
                XT.append(Xtemp)
                Ytemp = np.abs(ytemp-y).argmin()
                YT.append(Ytemp)
                LVTtemp = g16.variables['LVT'][Ytemp, Xtemp, j]
                LVTi.append(LVTtemp)
            LVT.append(LVTi)

        LVT = np.array(LVT)

        # Retrieve Relative humidity data
        g16_data_fileM = []
        xscanM = []
        yscanM = []
        LVM = []

        for i in range(0, len(nc_filesM)):
            g16_dataM = nc_filesM[i]
            g16_data_fileM.append(g16_dataM)
            g16M = Dataset(g16_data_fileM[i], 'r')
            xtempM = g16M.variables['x'][:]
            xscanM.append(xtempM)
            ytempM = g16M.variables['y'][:]
            yscanM.append(ytempM)

            XM = []
            YM = []
            LVMi = []
            for j in range(0, len(P)):
                XtempM = np.abs(xtempM-x).argmin()
                XM.append(XtempM)
                YtempM = np.abs(ytempM-y).argmin()
                YM.append(YtempM)
                LVMtemp = g16M.variables['LVM'][YtempM, XtempM, j]
                LVMi.append(LVMtemp)
            LVM.append(LVMi)

        LVM = np.array(LVM)

        # Change pressure units to Pa and Temperature to K
        P = 100*P
        LVT = LVT-273.15

        # Constants needed and integrand
        rho_w = 1000  # kg/m**3
        g = 9.81  # m/s**2
        C = (-1)/(rho_w*g)
        # ev = 100*6.11*LVM*10**((7.5*LVT)/(LVT+237.15)) # Partial water vapour pressure in Pa
        ev = 100*6.112*LVM*np.exp((17.67*LVT)/(LVT+243.5))
        q = (0.622*ev)/(P-0.378*ev)  # Specific humdity
        # Complete integrand multiplied by 1000 to get the PWV in mm.
        f = 1000*C*q

        # Numerical integration
        PWV = []
        for j in range(0, len(nc_filesT)):
            integral = 0
            for i in range(P_minb+1, P_maxb+1):
                integral = integral + (P[i]-P[i-1])*((f[j, i]+f[j, i-1])/2)
            PWV.append(integral)

        PWV = np.asarray(PWV)

        out_date = date

        if plot:
            # Plot and save data
            fig = plt.figure(figsize=(15, 10))
            ax = fig.add_subplot(111)
            ax.plot(hour, PWV, 'bo', ms=4)
            plt.title('Precipitable Water Vapor at zenith, {} on {}'.format(
                site, day[1]), fontsize=26)
            plt.xticks(rotation='vertical', fontsize=24)
            plt.yticks(fontsize=24)
            ax.set_xlabel("Date", color="C0", fontsize=24)
            ax.set_ylabel("PWV (mm)", color="C0", fontsize=24)
            every_nth = 6
            for n, label in enumerate(ax.xaxis.get_ticklabels()):
                if n % every_nth != 0:
                    label.set_visible(False)
            for n, label in enumerate(ax.xaxis.get_ticklines()):
                if n % every_nth != 0:
                    label.set_visible(False)
            plt.tight_layout()
            plt.show()
            fig.savefig('PWV_at_zenith_{}_{}.png'.format(site, day[1]))

        if csv:
            np.savetxt('PWV_zenith_{}_{}.csv'.format(site, day[1]), np.column_stack((date, PWV)),
                       delimiter=',', fmt='%s', header='Time,PWV', comments='')

    return

def los_pwv_data_generator(short_name, ut_dates, centroided_sources, interp_type='linear', force_output_path=''):
    """Runs fyodor_pwv and interpolates PWV measurements onto same time grid as observations for a chosen target.

    :param short_name: short name of the target
    :type short_name: str
    :param ut_dates: list of UT dates on which to generate PWV data. You must have a directory in your Calibrations/PWV/ directory for each UT data containing the LVMP and LVTP data from NOAA CLASS. 
    :type ut_dates: list of str
    :param centroided_sources: dataframe of source centroids
    :type centroided_sources: pandas DataFrame
    :param interp_type: 'linear' or 'cubicspline', the interpolation method used to interpolate from NOAA CLASS time stamps to Mimir time stamps, defaults to 'linear'
    :type interp_type: str, optional
    :param force_output_path: user-chosen directory to use in place of the default ~/Documents/PINES_analysis_toolkit directory for analysis, defaults to ''
    :type force_output_path: str, optional
    """
    
    # Get the user's pines_analysis_toolkit path
    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()

    sources = get_source_names(centroided_sources)
    #Find the coordinates of the target using the PINES sample excel file.
    pines_sample_path = pines_path/('Misc/PINES sample.xlsx')
    sample_df = pd.read_excel(pines_sample_path)
    boo_arr = [sources[0].split(' ')[1] in str(sample_df['Short Name'][i]) for i in range(len(sample_df))]
    target_row = np.where(boo_arr)[0][0]
    #target_row = np.where(np.array(sample_df['Short Name']) == sources[0])[0][0]
    target_ra = sample_df['RA (deg)'][target_row]
    target_dec = sample_df['Dec (deg)'][target_row]

    #Loop over dates and generate Fyodor data.
    for date in ut_dates:
        fyodor_data_path = pines_path/('Calibrations/PWV/'+date)
        fyodor_pwv(short_name, fyodor_data_path, P_min=775, P_max=300, RA=target_ra, Dec=target_dec, csv=True, force_output_path=force_output_path)
    
    #Now interpolate PWV data onto grid of PINES times.
    #Start by grabbing times and pwv values from fyodor output.
    time_strs = []
    pwv = []
    for date in ut_dates:
        fyodor_data_path = pines_path/('Objects/'+short_name+'/pwv/')
        #Read in the output csv from the previous step
        csv_path = fyodor_data_path/('PWV_los_'+date+'.csv')
        df = pd.read_csv(csv_path)
        time_strs.extend(df['Time'])
        pwv.extend(df['PWV'])

    #Convert the times to BJD TDB. 
    time_strs = np.array(time_strs)
    fyodor_dates = np.array(
        [datetime.strptime(i, '%Y-%m-%d %H:%M:%S') for i in time_strs])
    fyodor_times_jd_utc = np.array([julian.to_jd(i) for i in fyodor_dates])  # Convert to JD UTC
    
    #Convert JD UTC to BJD TDB
    fyodor_times_bjd_tdb = jd_utc_to_bjd_tdb(fyodor_times_jd_utc, target_ra, target_dec)

    pwv = np.array(pwv)

    #Remove any NaN entries
    fyodor_times = fyodor_times_bjd_tdb[~np.isnan(pwv)]
    pwv = pwv[~np.isnan(pwv)]

    #Get the times of the PINES observations.
    centroided_sources.columns = centroided_sources.columns.str.strip()
    pines_times = np.array(centroided_sources['Time BJD TDB'])

    # Interpolate Fyodor full disk pwv data onto the PINES data grid.
    if interp_type == 'linear':
        breakpoint()
        f1 = interpolate.interp1d(fyodor_times, pwv, kind='nearest')
        pwv_interp = f1(pines_times[1:])

    elif interp_type == 'cubicspline':
        tck = interpolate.splrep(fyodor_times, pwv)
        pwv_interp = interpolate.splev(pines_times, tck, der=0)

    #Write out interpolated PWV data to the object's pwv/ directory. 
    pwv_output_dict = {'Time BJD TDB': pines_times, 'PWV': pwv_interp}
    pwv_output_df = pd.DataFrame(pwv_output_dict)
    pwv_output_path = pines_path / \
        ('Objects/'+short_name+'/pwv/'+short_name+'_fyodor_pwv.csv')
    pwv_output_df.to_csv(pwv_output_path, index=0)
    
def readBT(file, R=None, npix=3.0, samp=2.5E7, waverange=None, air=False, DF=0.0,
           verbose=False, wave=None, flam=None):
    """
    Returns wavelength in microns and flux in erg s^-1 cm^-2 A^-1
    
    Arguments
    file - name of file to read from, ignored if wave and flam are set
    R    - resolution to convolve spectrum to, default to not change resolution
    npix - number of pixels per resolution element of output spectrum, most spectrographs use 2-4
    samp - sampling to interpolate wavelength grid to for convolution, must be higher than original sampling
    waverange - two-element array with start and end wavelengths of output ***IN MICRONS***
    air  - set to True to convert output wavelengths from vacuum to air
    DF   - logarithmic offset, 0.0 for Public Veyette models, -8.0 for other raw PHOENIX output
    verbose - set to True to print timing info
    wave - input wavelength array to process instead of reading in from file
    flam - input flux array to process instead of reading from file
    """
    
    #Default to no broadening
    if R is None:
        asIs = True
    else:
        asIs = False
    
    if waverange is None and not asIs:
        raise ValueError("You really don't want to run this without setting a waverange or reading in asIs. "
              "Your computer will probably freeze.")
    
    if verbose:
        stime = time.time()
        print('Starting File: ' + str(file))
    
    #NOT IMPLEMENTED YET
    # if vsini == 0:
        # vsini = 3.0 # km/s
        # rotBroad = 0
    # else:
        # rotbroad = 1
    
    #Read in data unless passed
    if verbose: stime1 = time.time()
    if wave is None or flam is None:
    
        readin = []
        if os.path.splitext(file)[1] == '.xz':
            import lzma
            openit = lzma.open
        else:
            openit = open
        with openit(file, 'rt') as file1:
            for line in file1:
                if '#' not in line:
                    subbed = re.sub(r'(\d)-(\d)', r'\1 -\2', line).replace('D','E')
                    readin.append(subbed.split())

        wave = []
        flam = []    
        for line in readin:
            if waverange[0] < float(line[0])/1e4 < waverange[1]:
                wave.append(float(line[0]))
                flam.append(float(line[1]))
        
        wave = np.array(wave)
        #flam = 10.0**(np.array(flam) + DF)
        flam = np.array(flam)
        wave, uwave = np.unique(wave, return_index=True)
        flam = flam[uwave]

    #breakpoint()
    if verbose: print('Reading file took ' + str(int(round(time.time() - stime1))) + ' seconds')
    
    #Return without processing if asIs is True
    if asIs:
        if waverange is not None:
            keep = (waverange[0] <= wave/1e4) & (wave/1e4 <= waverange[1])
            wave = wave[keep]
            flam = flam[keep]
    
        if verbose:
            etime = time.time()
            print('Finished in ' + str(int(round(etime-stime))) + ' seconds.')
        return wave, flam
        
    ss = 10.0 * npix
    
    #Default to full data range
    if waverange is None: waverange = np.array([wave[ss+1], wave[-1*ss - 2]]) * 1e-4
    
    #Interpolate to higher sampling
    if verbose: stime1 = time.time()
    sampBump = math.ceil(samp / R / npix)
    iWaveStart = waverange[0] * 1e4
    iWaveEnd = waverange[1] * 1e4
    #if np.min(wave) > iWaveStart or np.max(wave) < iWaveEnd:
    #    print('Oh no! Not enough wavelength range. File: ' + file)
    iWaveFull = iWaveStart * np.exp( np.arange( samp * np.log(iWaveEnd/iWaveStart) ) / samp )
    w = (iWaveFull > np.min(wave)) & (iWaveFull < np.max(wave))
    iWaveOffset = np.where(w)[0][0]
    iWave = iWaveFull[w]
    ww = (wave >= 0.9*iWaveStart) & (wave <= 1.1*iWaveEnd)
    tck = interpolate.splrep(wave[ww], flam[ww], s=0)
    iflam = np.array(interpolate.splev(iWave, tck, der=0))


    if verbose: print('Interpolation took ' + str(int(round(time.time() - stime1))) + ' seconds')
    
    #Convolve with spectrograph resolution
    if verbose: stime1 = time.time()
    fwhm = samp / R
    stdv = fwhm / ( 2.0 * math.sqrt( 2.0 * math.log(2.0) ) )
    #Create Kernel and normalize
    lsf  = np.array(scipy.signal.gaussian(int(4.0 * fwhm), stdv))
    lsf /= lsf.sum()
    #Convolve
    specflam = scipy.signal.fftconvolve(iflam, lsf, mode='same')
    if verbose: print('Convolution took ' + str(int(round(time.time() - stime1))) + ' seconds')
    
    #Integrate over sampBump to get npix per resolution element
    if verbose: stime1 = time.time()
    i1Full  = ( float(sampBump) * np.arange( np.floor(
              (np.size(iWaveFull)-1.0) / sampBump ) ) )        # beginning index of each integral
    i1sampFull = i1Full + sampBump
    w = (iWaveFull[i1Full.astype('int').tolist()] > np.min(wave)) & (iWaveFull[i1sampFull.astype('int').tolist()] < np.max(wave))
    i1     = (i1Full.astype('int')[w] - iWaveOffset).tolist()
    i1samp = (i1sampFull.astype('int')[w] - iWaveOffset).tolist()
    trapA   = ( ( iWave[1:] - iWave[0:-1] ) *
              0.5 * ( specflam[0:-1] + specflam[1:] ) )         # trapezoidal area under each samp point
    intFlux = np.sum((trapA[int(i1[0]):int(i1[-1]+sampBump)]
                      ).reshape(len(i1), sampBump), axis=1)     # sum over sampBump points for each pixel
    #import pdb ; pdb.set_trace()
    delLam  = iWave[i1samp] - iWave[i1]   # delta lambda
    intWave = ( (iWave[i1samp]**2  - iWave[i1]**2)
              / ( 2.0 * delLam ) )                              # mean wavelength of each pixel
    intflam = intFlux / delLam
    nSpec   = len(intWave)
    if verbose: print('Integration took ' + str(int(round(time.time() - stime1))) + ' seconds')
    
    if air: intWave /= (1.0 + 2.735182e-4 + 131.4182 / intWave**2 + 2.76249e8 / intWave**4)
    
    if verbose:
        etime = time.time()
        print('Finished in ' + str(int(round(etime-stime))) + ' seconds.')
    
    return intWave/1e4, intflam

def pwv_spectrum_time_series(short_name, reference_id, band='J', force_output_path=''):
    plate_scale = 0.579

    # Get the user's pines_analysis_toolkit path
    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()
    
    if band == 'J':
        waverange = [1.125, 1.375]
        waverange_readBT = [0.5, 2.0]
        filter_path = pines_path/('Misc/'+'mimir_'+band+'_bandpass.txt')
        filter_df = pd.read_csv(filter_path, sep='\t', names=['Wavelength', 'Throughput'])
        mimir_w =  np.array(filter_df['Wavelength'])
        mimir_thru = np.array(filter_df['Throughput'])/100
        inds = np.where((waverange[0]-0.01 < mimir_w) & (mimir_w < waverange[1]+0.01))[0]
        mimir_w = mimir_w[inds]
        mimir_thru = mimir_thru[inds]
    else:
        print('Have not implemented anything other than J band, returning!')
        return
    
    object_path = pines_path/('Objects/'+short_name)
    centroid_df = pines_log_reader(object_path/('sources/target_and_references_centroids.csv'))
    airmass = np.array(centroid_df['Airmass'])
    seeing = np.array(centroid_df['Seeing'])
    zenith_angle = np.arccos(1/airmass)*180/np.pi
    bt_settl_path = pines_path/('Calibrations/BT-Settl')

    atran_path = pines_path/('Calibrations/PWV/ATRAN_spectra')
    #Get info about the reference star
    source_df = pd.read_csv(object_path/('sources/target_and_references_source_detection.csv'))
    entry = source_df.iloc[reference_id]
    temp = entry['Teff']

    #Get the reference star's spectrum. Units are erg cm^-2 s^-1 A^-1
    if int(temp/100) <= 25:
        bt_settl_filename = 'lte0'+str(int(temp/100))+'-4.5-0.0.BT-Settl.7.dat.txt'
    else:
        bt_settl_filename = 'lte0'+str(int(temp/100))+'-4.5-0.0a+0.0.BT-NextGen.7.dat.txt'

    print('Source {} temp. = {} K, using {} BT-Settl model.'.format(reference_id, temp, bt_settl_filename))
    spectrum_w, spectrum_flam = readBT(bt_settl_path/bt_settl_filename, R=100, waverange=waverange_readBT)

    #Convert from ergs, [J cm^-2 s^-1 A^-1]
    spectrum_flam = spectrum_flam / 1e7

    #Multiply by wavelength in Angstrom, [J cm^-2 s^-1]
    spectrum_flam = spectrum_flam * spectrum_w * 1e4 / (100*3)

    #Multiply by the area of Perkins in cm^2, [J s^-1]
    A = (1-0.05)*np.pi*(1.83*100/2)**2
    spectrum_flam = spectrum_flam * A

    #Multiply by exposure time, [J]
    exptime = 20 
    spectrum_flam = spectrum_flam * exptime

    #Divide by Joules/photon, [photons]
    h = 6.62607015e-34
    c = 2.99792458e8
    E_phot = h*c/(spectrum_w*1e-6) #[J/photon]
    spectrum_flam = spectrum_flam/E_phot

    # #Multiply by (R/D)**2
    # if reference_id == 1:
    #     R_star = 0.8 #solar radii
    #     p = 0.665 #mas
    # elif reference_id == 5:
    #     R_star = 0.3 #solar radii 
    #     p = 5.269 #mas
    # else:
    #     breakpoint()
    # D = 1/(p/1000) * 3.086e16
    # R = R_star * 7e8
    # print((R/D)**2)
    # spectrum_flam = spectrum_flam *  (R/D)**2 
    
    #Get the PWV data.
    pines_pwv_path = object_path/('pwv/'+short_name+'_fyodor_pwv.csv')
    pines_pwv_df = pd.read_csv(pines_pwv_path)
    pines_times = np.array(pines_pwv_df['Time BJD TDB'])
    pines_pwv = np.array(pines_pwv_df['PWV'])

    integrated_response = np.zeros(len(pines_times))
    seeing_factors = np.zeros(len(pines_times))
    pwv_factors = np.zeros(len(pines_times))
    for i in range(len(pines_times)):
        #Round zenith angle to the nearest 5 degrees. 
        z = round(zenith_angle[i]/5)*5
        #Round PWV to the nearest 0.5
        pwv = round(pines_pwv[i]*2)/2

        #Get the appropriate ATRAN sky spectrum given z and pwv
        atran_file = 'atran_{:1.2f}_{:1d}.dat'.format(pwv, z)
        #atran_file = 'atran_{:1.2f}_35.dat'.format(pwv)

        atran_df = pd.read_csv(atran_path/atran_file, delim_whitespace=True, names=['Wavelength', 'Transmission'])
        atran_w = np.array(atran_df['Wavelength'])
        atran_tran = np.array(atran_df['Transmission'])
        
        #Interpolate Mimir throughput, sky transmission, and stellar spectrum onto the same wavelength grid. 
        interp_w = np.linspace(waverange[0], waverange[1], 1000)

        f_mimir = scipy.interpolate.interp1d(mimir_w, mimir_thru)
        f_sky = scipy.interpolate.interp1d(atran_w, atran_tran)
        f_star = scipy.interpolate.interp1d(spectrum_w, spectrum_flam)

        interp_mimir = f_mimir(interp_w)
        interp_sky = f_sky(interp_w)
        interp_star = f_star(interp_w)

        seeing_sigma = seeing[i]/(plate_scale)
        seeing_factors[i] = norm(0, seeing_sigma).cdf(4)-norm(0, seeing_sigma).cdf(-4)

        #Sum over wavelength and correct for seeing effects
        pwv_factors[i] = np.sum(interp_star*interp_sky*interp_mimir)
        integrated_response[i] = pwv_factors[i]
        #integrated_response[i] = pwv_factors[i] * seeing_factors[i] * 0.058 #Extra 0.058 to account for Mimir throughput
    
    integrated_response = integrated_response/np.mean(integrated_response)
    integrated_response = convolve(integrated_response, Gaussian1DKernel(31), boundary='extend')

    return integrated_response
