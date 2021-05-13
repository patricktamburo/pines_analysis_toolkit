from netCDF4 import Dataset                                             
import os                                                               
import numpy as np                                                    
import time                                                             
import matplotlib.pyplot as plt  
import matplotlib.patches as mpatches
import glob
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import EarthLocation, AltAz
from astropy.time import Time
from astropy.utils.data import clear_download_cache
from astropy.constants import R_earth
import requests
from bs4 import BeautifulSoup
import wget
from multiprocessing.pool import ThreadPool
import pdb 
from pines_analysis_toolkit.utils import pines_dir_check, short_name_creator

"""
    THIS IS AN ADAPTATION OF THE FYODOR PACKAGE TO WORK WITH PINES DATA https://github.com/Erikmeier18/fyodor
"""

__all__ = ['rename_nc', 'download_nc', 'pwv']

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
    os.chdir(str(directory))
    full_direc = os.listdir()
    for i in range(len(full_direc)):
        try: 
            tmp = int(full_direc[i].split('.')[0])
            dst = full_direc[i].split('.')[1] + '.' + full_direc[i].split('.')[2]
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

def download_nc(url, directory, date, n_threads=5):
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

    if len(str(date[7])) == 1:
        subs = '{}'.format(date[0]) + '00' + '{}'.format(date[7])
    elif len(str(date[7])) == 2:
        subs = '{}'.format(date[0]) + '0' + '{}'.format(date[7])
    elif len(str(date[7])) == 3:
        subs = '{}'.format(date[0]) + '{}'.format(date[7])

    #Modified this so that we know how many files to download before we start.
    FileName = [i for i in Content if subs in i]
    FileName = [path for path in FileName if not os.path.exists('./'+ path)]

    urls = []
    for i in range(0, len(FileName)):
        urlstemp = url + "/" + FileName[i] 
        urls.append(urlstemp)

    print('There are %s files to download' % len(urls))

    for ii in range(0, len(urls), n_threads):
        slice_item = slice(ii, ii + n_threads, 1)

        files = ThreadPool(len(urls[slice_item])).imap_unordered(download_file_wget, urls[slice_item])
        # Need to do this print so that it waits before starting the next batch.
        for f in files:
            print('%s complete' % f)

    print("All files downloaded!")
    

def pwv(target, directory, P_min, P_max, line_of_sight='target', RA=None, Dec=None, plot=False, csv=False):
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

    # Access working directory
    os.chdir(directory)
    #nc_filesT = glob.glob('*OR_ABI-L2-LVTPF*') #MODIFIED TO WORK WITH CONUS DATA
    nc_filesT = glob.glob('*OR_ABI-L2-LVTP*') 

    nc_filesT = sorted(nc_filesT)
    #nc_filesM = glob.glob('*OR_ABI-L2-LVMPF*')
    nc_filesM = glob.glob('*OR_ABI-L2-LVMP*')
    nc_filesM = sorted(nc_filesM) 

    # Open the first file to retrieve earth parameters
    Proj_info = Dataset(nc_filesT[0],'r')
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
        epochtemp = 946728000+ int(t[i])
        epoch.append(epochtemp)
        datetemp = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(epoch[i]))
        date.append(datetemp)
        daytemp = time.strftime("%d-%m-%Y", time.gmtime(epoch[i]))
        day.append(daytemp)
        hourtemp = time.strftime("%H:%M:%S", time.gmtime(epoch[i]))
        hour.append(hourtemp)
    
    #Check if the output already exists, if so return. 
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    out_date = day[1][-4:]+day[1][3:5]+day[1][0:2]
    #Make the object's pwv directory if it doesn't already exist.
    if not os.path.isdir(pines_path/('Objects/'+short_name+'/pwv/')):
        os.mkdir(pines_path/('Objects/'+short_name+'/pwv/'))
    output_path = pines_path/('Objects/'+short_name+'/pwv/PWV_los_{}.csv'.format(out_date))
    if os.path.exists(output_path):
        print('PWV output already exists for {}, returning.'.format(out_date))
        return 

    # Use astropy.time to keep format for target coordinates:
    times = Time(date, format='iso', scale='utc')

    # Barometric formula
    p0 = P[0] #hPa
    R_D = 287 #Jkg-1K-1
    g = 9.81 #m/s2
    T_s = 288 #K
    L = -0.0065 #K/m

    h = (T_s/L)*((P/p0)**(-L*R_D/g)-1)*u.m

    e = np.sqrt((r_eq**2-r_pol**2)/(r_eq**2)) #Eccentricity

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
        INDEX = np.ravel(np.where(Aa.alt.degree<0)) #ORIGINAL USED VALUES WERE WHERE ALT WAS >30. CHANGED FOR PERKINS.
        INDEXP = np.ravel(np.where(Aa.alt.degree>0))
    
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

        for i in range(0,len(Alt)):
            delta_xi = []
            for j in h:
                delta_xt = j/np.tan(Alt[i])
                delta_xt = delta_xt*u.m**-1
                delta_xi.append(delta_xt)
            delta_x.append(delta_xi)
        delta_x = np.array(delta_x)

        for i in range(0,len(Az)):
            d_latt = delta_x[i,]*np.cos(Az[i])
            d_lont = delta_x[i,]*np.sin(Az[i])
            d_lat.append(d_latt)
            d_lon.append(d_lont)
        d_lat = np.array(d_lat)
        d_lon = np.array(d_lon)

        # Compute latitude and longitude of projection points
        lat_proj = []
        lon_proj = []
        for i in range(0, len(Alt)):
            obs_latt = loc.lat.degree + raddeg*(np.arctan(d_lat[i,]/R_earth*u.m**1)*u.rad**-1) 
            obs_lont = loc.lon.degree + raddeg*(np.arctan(d_lon[i,]/R_earth*u.m**1)*u.rad**-1)
            lat_proj.append(obs_latt)
            lon_proj.append(obs_lont)
        lat_proj = np.array(lat_proj)
        lon_proj = np.array(lon_proj)
        rad = (np.pi)/180

        lat_proj_rad = rad*lat_proj
        lon_proj_rad = rad*lon_proj
        lambda_0 = rad*lon_origin
    
        #T ransform into scan angles
        lat_origin = np.arctan(((r_pol**2)/(r_eq**2))*np.tan(lat_proj_rad))
    
        r_c = r_pol/(np.sqrt(1-(e**2)*(np.cos(lat_origin))**2))

        s_x = H -r_c*np.cos(lat_origin)*np.cos(lon_proj_rad-lambda_0)
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
                Xtemp = np.abs( xtemp-x[i,j]).argmin()
                Xi.append(Xtemp)
                Ytemp = np.abs( ytemp-y[i,j]).argmin()
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
                XtempM = np.abs( xtempM-x[i,j]).argmin()
                Xi.append(XtempM)
                YtempM = np.abs( ytempM-y[i,j]).argmin()
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
        g = 9.81  #m/s**2
        C = (-1)/(rho_w*g)
        ev = 100*6.112*LVM*np.exp((17.67*LVT)/(LVT+243.5))  # Partial water vapour pressure in Pa
        q = (0.622*ev)/(Pi-0.378*ev)  #Specific humdity
        f = 1000*C*q  #Complete integrand multiplied by 1000 to get the PWV in mm.
    
        # Numerical integration
        PWV = [] 
        for j in range(0, len(LVT)):
            integral = 0
            for i in range(1, len(Pi)): 
                integral = integral +(Pi[i]-Pi[i-1])*((f[j,i]+f[j,i-1])/2)  
            PWV.append(integral)

        PWV = np.asarray(PWV)
    
        if plot:
            # Plot and save data
            fig = plt.figure(figsize=(15,10))
            ax=fig.add_subplot(111)
            ax.plot(HOUR, PWV, 'bo', ms=4)
            plt.title('PWV along line of sight, {} on {}'.format(site, day[1]), fontsize=26)
            plt.xticks(rotation='vertical', fontsize=24)
            plt.yticks(fontsize=24)
            ax.set_xlabel("Date", color="C0", fontsize =24)
            ax.set_ylabel("PWV (mm)", color="C0", fontsize =24)
            RA_patch = mpatches.Patch(color='white', label='RA: {} degrees'.format(RA))
            Dec_patch = mpatches.Patch(color='white', label='Dec: {} degrees'.format(Dec))
            every_nth = 6
            for n, label in enumerate(ax.xaxis.get_ticklabels()):
                if n % every_nth != 0:
                    label.set_visible(False)
            for n, label in enumerate(ax.xaxis.get_ticklines()):
                if n % every_nth != 0:
                    label.set_visible(False)
            plt.tight_layout()
            plt.legend(handles=[RA_patch, Dec_patch], loc='lower right', fontsize=22)
            plt.show()
            fig.savefig('PWV_line_of_sight_{}_{}.png'.format(site, day[1]))
        
        if csv:
            np.savetxt(output_path, np.column_stack((date, PWV)), 
                   delimiter=',' , fmt = '%s', header= 'Time,PWV', comments='')

    # Computes PWV at zenith
    elif line_of_sight == 'zenith':
    
        # Transform latitude and longitude into scan angles
        rad = np.pi/180
        lambda_0 = rad*lon_origin
        obs_lat_rad = rad*latt
        obs_lon_rad = rad*lont
    
        lat_origin = np.arctan(((r_pol**2)/(r_eq**2))*np.tan(obs_lat_rad))

        r_c = r_pol/(np.sqrt(1-(e**2)*(np.cos(lat_origin))**2))

        s_x = H -r_c*np.cos(lat_origin)*np.cos(obs_lon_rad-lambda_0)
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
                Xtemp = np.abs( xtemp-x).argmin() 
                XT.append(Xtemp)
                Ytemp = np.abs( ytemp-y).argmin()
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
                XtempM = np.abs( xtempM-x).argmin()
                XM.append(XtempM)
                YtempM = np.abs( ytempM-y).argmin()
                YM.append(YtempM)
                LVMtemp = g16M.variables['LVM'][YtempM, XtempM, j]
                LVMi.append(LVMtemp)
            LVM.append(LVMi)

        LVM = np.array(LVM)

        # Change pressure units to Pa and Temperature to K
        P = 100*P 
        LVT = LVT-273.15
    
        # Constants needed and integrand
        rho_w = 1000 # kg/m**3
        g = 9.81 #m/s**2
        C = (-1)/(rho_w*g)
        #ev = 100*6.11*LVM*10**((7.5*LVT)/(LVT+237.15)) # Partial water vapour pressure in Pa
        ev = 100*6.112*LVM*np.exp((17.67*LVT)/(LVT+243.5))
        q = (0.622*ev)/(P-0.378*ev) #Specific humdity
        f = 1000*C*q #Complete integrand multiplied by 1000 to get the PWV in mm.

        # Numerical integration
        PWV = [] 
        for j in range(0, len(nc_filesT)): 
            integral = 0
            for i in range(P_minb+1, P_maxb+1):
                integral = integral +(P[i]-P[i-1])*((f[j,i]+f[j,i-1])/2)  
            PWV.append(integral)

        PWV = np.asarray(PWV)

        out_date=date
    
        if plot:
            # Plot and save data
            fig = plt.figure(figsize=(15,10))
            ax=fig.add_subplot(111)
            ax.plot(hour, PWV, 'bo', ms=4)
            plt.title('Precipitable Water Vapor at zenith, {} on {}'.format(site, day[1]), fontsize=26)
            plt.xticks(rotation='vertical', fontsize=24)
            plt.yticks(fontsize=24)
            ax.set_xlabel("Date", color="C0", fontsize =24)
            ax.set_ylabel("PWV (mm)", color="C0", fontsize =24)
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
            np.savetxt('PWV_zenith_{}_{}.csv'.format(site,day[1]), np.column_stack((date, PWV)), 
                   delimiter=',' , fmt = '%s', header= 'Time,PWV', comments='')
    
    return 
