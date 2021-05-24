#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 21:57:40 2021

@authors: davidgracia, Patrick Tamburo 
"""

import requests
import json
from pines_analysis_toolkit.astrometry.astrometry_client import Client
import time
import pdb 

astrometry = Client()
    
def upload_file(apikey, filename, header):
    
    """ uploads multiple file to astrometry.net
    
        number = number of files to upload
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
            
