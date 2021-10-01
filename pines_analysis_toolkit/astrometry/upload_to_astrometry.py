import requests
from pines_analysis_toolkit.astrometry.astrometry_client import Client
import time

astrometry = Client()
    
def upload_file(apikey, filename, header):
    """Uploads files to astrometry.net. 

    :param apikey: api key of your account on astrometry.net
    :type apikey: str
    :param filename: path to the file to upload
    :type filename: pathlib.PosixPath
    :param header: header of the file you want to upload
    :type header: astropy.io fits header 
    """
    
    astrometry.login(apikey)
    
    tel_ra_str = header['TELRA'].split(':')
    tel_ra_deg = int(tel_ra_str[0])*15 + int(tel_ra_str[1])*(15/60) + float(tel_ra_str[2])*(15/3600)

    tel_dec_str = header['TELDEC'].split(':')
    
    if tel_dec_str[0][0] == '+':
        tel_dec_deg = int(tel_dec_str[0]) + int(tel_dec_str[1])/60 + float(tel_dec_str[2])/3600
    elif tel_dec_str[0][0] == '-':
        tel_dec_deg = int(tel_dec_str[0]) - int(tel_dec_str[1])/60 - float(tel_dec_str[2])/3600

    sub_id = astrometry.upload(filename, scale_units='arcsecperpix', scale_lower=0.575, scale_upper=0.585, center_ra=tel_ra_deg, center_dec=tel_dec_deg, radius=0.116)['subid']
    print('submission_id: ' , str(sub_id))
    
    while astrometry.sub_status(sub_id)['jobs'] == [] or astrometry.sub_status(sub_id)['jobs'] == [None]:
        time.sleep(1)
        
    job_id = astrometry.sub_status(sub_id)['jobs'][0]
    print('job_id: ', astrometry.sub_status(sub_id)['jobs'])
    
    while astrometry.job_status(job_id) == 'solving':
        time.sleep(1)

    outfile = filename.parent/(filename.name.split('.fits')[0]+'_new_image.fits')
    
    url = 'http://nova.astrometry.net/new_fits_file/' + str(job_id)
    
    r = requests.get(url, allow_redirects=True)
    open(outfile, 'wb').write(r.content)
            
