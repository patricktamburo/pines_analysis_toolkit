#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 21:57:40 2021

@author: davidgracia
"""

import requests
import json
from pines_analysis_toolkit.astrometry.astrometry_client import Client
import time
import pdb 

astrometry = Client()

def login():
    
    """ Login in to astrometry.net and get a session ID
    """
    
    apikey = 'mioghoxehfwdxyim'
    astrometry.login(apikey)
    

# def upload_file(apikey, filename):
    
#     """ uploads file to astrometry.net
#     """
    
#     #apikey = 'mioghoxehfwdxyim'
#     astrometry.login(apikey)
    
#     astrometry.upload(filename) 
    
def upload_file(apikey, filename):
    
    """ uploads multiple file to astrometry.net
    
        number = number of files to upload
    """
    
    astrometry.login(apikey)
    
    # filename = '/Users/davidgracia/Downloads/Astronomy Research/20201105.100_red_no_nan.fits'
    # beg = filename[:57]
    # end = filename[60:]
    # num = 100
    
        
    
    sub_id = astrometry.upload(filename)['subid']
    print(' ')
    print('sub_id: ' , str(sub_id))
    print(' ')
    
    while astrometry.sub_status(sub_id)['jobs'] == [] or astrometry.sub_status(sub_id)['jobs'] == [None]:
        time.sleep(5)
        
    job_id = astrometry.sub_status(sub_id)['jobs'][0]
    print(' ')
    print('job_id: ', astrometry.sub_status(sub_id)['jobs'])
    print(' ')
    
    while astrometry.job_status(job_id) == 'solving':
        time.sleep(1)
    
    outfile = filename.split('.fits')[0] + '_new_image.fits'
    
    url = 'http://nova.astrometry.net/new_fits_file/' + str(job_id)
    
    r = requests.get(url, allow_redirects=True)
    open(outfile, 'wb').write(r.content)
            
        

def download_file(version, job_id, filename):
    
    """ Downoald requested file from astrometry
    
        for wcs, version = 'wcs'
    
        for new-image, version = 'new_fits'
        
        for corr, version = 'corr_file'
    """
    
    url = "https://nova.astrometry.net/" + version + "_file/" + job_id
    
    r = requests.get(url, allow_redirects=True)
    open(filename, 'wb').write(r.content)



#functions to use in the command line to check status on job/submission
"""
astrometry.myjobs()
astrometry.job_status(job_id)
astrometry.sub_status(sub_id)
astrometry.submission_images(subid)
"""
