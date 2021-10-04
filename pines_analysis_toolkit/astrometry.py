from __future__ import print_function

import pines_analysis_toolkit as pat 

from astropy.convolution import interpolate_replace_nans, Gaussian2DKernel
from astropy.io import fits 
from astropy.wcs import WCS

from glob import glob 
from natsort import natsorted 
import numpy as np
import os 
import time
from pathlib import Path
import pandas as pd 
import requests

import sys
import base64

from urllib.parse import urlencode, quote
from urllib.request import urlopen, Request
from urllib.error import HTTPError

from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.application  import MIMEApplication
from email.encoders import encode_noop

import json

def json2python(data):
    try:
        return json.loads(data)
    except:
        pass
    return None
python2json = json.dumps

class MalformedResponse(Exception):
    pass

class RequestError(Exception):
    pass

class Client(object):
    default_url = 'http://nova.astrometry.net/api/'

    def __init__(self,
                 apiurl = default_url):
        self.session = None
        self.apiurl = apiurl

    def get_url(self, service):
        return self.apiurl + service

    def send_request(self, service, args={}, file_args=None):
        '''
        service: string
        args: dict
        '''
        if self.session is not None:
            args.update({ 'session' : self.session })
        #print('Python:', args)
        json = python2json(args)
        #print('Sending json:', json)
        url = self.get_url(service)
        #print('Sending to URL:', url)

        # If we're sending a file, format a multipart/form-data
        if file_args is not None:
            import random
            boundary_key = ''.join([random.choice('0123456789') for i in range(19)])
            boundary = '===============%s==' % boundary_key
            headers = {'Content-Type':
                       'multipart/form-data; boundary="%s"' % boundary}
            data_pre = (
                '--' + boundary + '\n' +
                'Content-Type: text/plain\r\n' +
                'MIME-Version: 1.0\r\n' +
                'Content-disposition: form-data; name="request-json"\r\n' +
                '\r\n' +
                json + '\n' +
                '--' + boundary + '\n' +
                'Content-Type: application/octet-stream\r\n' +
                'MIME-Version: 1.0\r\n' +
                'Content-disposition: form-data; name="file"; filename="%s"' % file_args[0] +
                '\r\n' + '\r\n')
            data_post = (
                '\n' + '--' + boundary + '--\n')
            data = data_pre.encode() + file_args[1] + data_post.encode()

        else:
            # Else send x-www-form-encoded
            data = {'request-json': json}
            #print('Sending form data:', data)
            data = urlencode(data)
            data = data.encode('utf-8')
            #print('Sending data:', data)
            headers = {}

        request = Request(url=url, headers=headers, data=data)

        try:
            f = urlopen(request)
            txt = f.read()
            #print('Got json:', txt)
            result = json2python(txt)
            #print('Got result:', result)
            stat = result.get('status')
            #print('Got status:', stat)
            if stat == 'error':
                errstr = result.get('errormessage', '(none)')
                raise RequestError('server error message: ' + errstr)
            return result
        except HTTPError as e:
            print('HTTPError', e)
            txt = e.read()
            open('err.html', 'wb').write(txt)
            print('Wrote error text to err.html')

    def login(self, apikey):
        args = { 'apikey' : apikey }
        result = self.send_request('login', args)
        sess = result.get('session')
        #print('Got session:', sess)
        if not sess:
            raise RequestError('no session in result')
        self.session = sess

    def _get_upload_args(self, **kwargs):
        args = {}
        for key,default,typ in [('allow_commercial_use', 'd', str),
                                ('allow_modifications', 'd', str),
                                ('publicly_visible', 'n', str),
                                ('scale_units', None, str),
                                ('scale_type', None, str),
                                ('scale_lower', None, float),
                                ('scale_upper', None, float),
                                ('scale_est', None, float),
                                ('scale_err', None, float),
                                ('center_ra', None, float),
                                ('center_dec', None, float),
                                ('parity',None,int),
                                ('radius', None, float),
                                ('downsample_factor', None, int),
                                ('positional_error', None, float),
                                ('tweak_order', None, int),
                                ('crpix_center', None, bool),
                                ('invert', None, bool),
                                ('image_width', None, int),
                                ('image_height', None, int),
                                ('x', None, list),
                                ('y', None, list),
                                ('album', None, str),
                                ]:
            if key in kwargs:
                val = kwargs.pop(key)
                val = typ(val)
                args.update({key: val})
            elif default is not None:
                args.update({key: default})
        
        #print('Upload args:', args)
        return args

    def url_upload(self, url, **kwargs):
        args = dict(url=url)
        args.update(self._get_upload_args(**kwargs))
        result = self.send_request('url_upload', args)
        return result

    def upload(self, fn=None, **kwargs):
        args = self._get_upload_args(**kwargs)
        file_args = None
        if fn is not None:
            try:
                f = open(fn, 'rb')
                file_args = (fn, f.read())
            except IOError:
                print('File %s does not exist' % fn)
                raise
        return self.send_request('upload', args, file_args)

    def submission_images(self, subid):
        result = self.send_request('submission_images', {'subid':subid})
        return result.get('image_ids')

    def overlay_plot(self, service, outfn, wcsfn, wcsext=0):
        from astrometry.util import util as anutil
        wcs = anutil.Tan(wcsfn, wcsext)
        params = dict(crval1 = wcs.crval[0], crval2 = wcs.crval[1],
                      crpix1 = wcs.crpix[0], crpix2 = wcs.crpix[1],
                      cd11 = wcs.cd[0], cd12 = wcs.cd[1],
                      cd21 = wcs.cd[2], cd22 = wcs.cd[3],
                      imagew = wcs.imagew, imageh = wcs.imageh)
        result = self.send_request(service, {'wcs':params})
        #print('Result status:', result['status'])
        plotdata = result['plot']
        plotdata = base64.b64decode(plotdata)
        open(outfn, 'wb').write(plotdata)
        #print('Wrote', outfn)

    def sdss_plot(self, outfn, wcsfn, wcsext=0):
        return self.overlay_plot('sdss_image_for_wcs', outfn,
                                 wcsfn, wcsext)

    def galex_plot(self, outfn, wcsfn, wcsext=0):
        return self.overlay_plot('galex_image_for_wcs', outfn,
                                 wcsfn, wcsext)

    def myjobs(self):
        result = self.send_request('myjobs/')
        return result['jobs']

    def job_status(self, job_id, justdict=False):
        result = self.send_request('jobs/%s' % job_id)
        if justdict:
            return result
        stat = result.get('status')
        if stat == 'success':
            result = self.send_request('jobs/%s/calibration' % job_id)
            #print('Calibration:', result)
            result = self.send_request('jobs/%s/tags' % job_id)
            #print('Tags:', result)
            result = self.send_request('jobs/%s/machine_tags' % job_id)
            #print('Machine Tags:', result)
            result = self.send_request('jobs/%s/objects_in_field' % job_id)
            #print('Objects in field:', result)
            result = self.send_request('jobs/%s/annotations' % job_id)
            #print('Annotations:', result)
            result = self.send_request('jobs/%s/info' % job_id)
            #print('Calibration:', result)

        return stat

    def annotate_data(self,job_id):
        """
        :param job_id: id of job
        :return: return data for annotations
        """
        result = self.send_request('jobs/%s/annotations' % job_id)
        return result

    def sub_status(self, sub_id, justdict=False):
        result = self.send_request('submissions/%s' % sub_id)
        if justdict:
            return result
        """
        return result.get('status')
        """
        return result

    def jobs_by_tag(self, tag, exact):
        exact_option = 'exact=yes' if exact else ''
        result = self.send_request(
            'jobs_by_tag?query=%s&%s' % (quote(tag.strip()), exact_option),
            {},
        )
        return result
        
def pines_astrometry(target, api_key, download_data=False):
    """Uploads reduced images for a target to astrometry.net, downloads solution image, and updates the image header with the astrometry.net wcs.

    :param target: long name of the target
    :type target: str
    :param api_key: he api key of your account on astrometry.net
    :type api_key: str
    :param download_data: whether or not to first download reduced data for the target from the PINES server, defaults to False
    :type download_data: bool, optional
    """
    pines_path = pat.utils.pines_dir_check()
    short_name = pat.utils.short_name_creator(target)

    if download_data:
        sftp = pat.utils.pines_login()
        pat.data.get_reduced_science_files(sftp, target)

    reduced_path = pines_path/('Objects/'+short_name+'/reduced/')
    reduced_files = np.array(natsorted(glob(str(reduced_path)+'/*red.fits')))

    kernel = Gaussian2DKernel(x_stddev=0.25)

    times = []
    for i in range(0, len(reduced_files)):
        t1 = time.time()
        filename = Path(reduced_files[i])
        print(short_name+', '+filename.name+', image '+str(i+1)+' of '+str(len(reduced_files)))

        #If the header does not have a HISTORY keyword (which is added by astrometry.net), process it. 
        header = fits.open(filename)[0].header

        if 'HISTORY' not in header:
            #Read in the image data, interpolate NaNs, and save to a temporary fits file. 
            #Astrometry.net does not work with NaNs in images. 
            original_image = fits.open(filename)[0].data
            image = interpolate_replace_nans(original_image, kernel=kernel)
            temp_filename = filename.parent/(filename.name.split('.fits')[0]+'_temp.fits')
            hdu = fits.PrimaryHDU(image, header=header)
            hdu.writeto(temp_filename, overwrite=True)

            #Upload the image with interpolated NaNs to astrometry.net. It will solve and download to a new image. 
            pat.astrometry.upload_to_astrometry.upload_file(api_key, temp_filename, header)

            #Try to donwload the solved image and open it. 
            try:
                #Grab the header of the astrometry.net solution image, and the original image data. 
                astrometry_image_path = filename.parent/(temp_filename.name.split('.fits')[0]+'_new_image.fits')
                wcs_header = fits.open(astrometry_image_path)[0].header
                wcs_hdu = fits.PrimaryHDU(original_image, header=wcs_header)

                #Save the original image data with the new wcs header. 
                output_filename = filename
                wcs_hdu.writeto(output_filename, overwrite=True)

            #If the try clause didn't work, that's because the processing on astrometry.net failed. 
            except:
                print('Astrometry solution failed!')
                #Add HISTORY keyword to header so it doesn't get processed in the future. 
                header['HISTORY'] = 'Astrometric solution failed.'
                fits.writeto(filename, original_image, header, overwrite=True)

            #Delete temporary files. 
            os.remove(temp_filename)
            os.remove(astrometry_image_path)
            time.sleep(1)
            t2 = time.time()
            times.append(t2-t1)
            print('Average completion time: {:3.1f} s.'.format(np.mean(times)))
            print('')

        #If the header DOES have a HISTORY keyword, skip it, it has already been processed. 
        else:
            print('Astrometric processing already complete for {}, skipping.'.format(filename.name))
            print('')

def source_detect_astrometry(api_key, filename):
    """ Uploads a single reduced image for a target to astrometry.net, downloads solution image, and updates the image header with the astrometry.net wcs. 

    :param api_key: the api key of your account on astrometry.net
    :type api_key: str
    :param filename: path to the source_detect_image
    :type filename: pathlib.PosixPath
    :raises RuntimeError: if astrometry solution fails
    """
    kernel = Gaussian2DKernel(x_stddev=0.25)

    #If the header does not have a HISTORY keyword (which is added by astrometry.net), process it. 
    header = fits.open(filename)[0].header

    if 'HISTORY' not in header:
        print('Uploading {} to astrometry.net.'.format(filename.name))
        #Read in the image data, interpolate NaNs, and save to a temporary fits file. 
        #Astrometry.net does not work with NaNs in images. 
        original_image = fits.open(filename)[0].data
        image = interpolate_replace_nans(original_image, kernel=kernel)
        temp_filename = filename.parent/(filename.name.split('.fits')[0]+'_temp.fits')
        hdu = fits.PrimaryHDU(image, header=header)
        hdu.writeto(temp_filename, overwrite=True)

        #Upload the image with interpolated NaNs to astrometry.net. It will solve and download to a new image. 
        pat.astrometry.upload_file(api_key, temp_filename, header)

        #Try to donwload the solved image and open it. 
        try:
            #Grab the header of the astrometry.net solution image, and the original image data. 
            astrometry_image_path = filename.parent/(temp_filename.name.split('.fits')[0]+'_new_image.fits')
            wcs_header = fits.open(astrometry_image_path)[0].header
            wcs_hdu = fits.PrimaryHDU(original_image, header=wcs_header)

            #Save the original image data with the new wcs header. 
            output_filename = filename
            wcs_hdu.writeto(output_filename, overwrite=True)

        #If the try clause didn't work, that's because the processing on astrometry.net failed. 
        except:
            raise RuntimeError('Astrometry solution failed! Use a different source detect image for choosing reference stars.')            

        #Delete temporary files. 
        os.remove(temp_filename)
        os.remove(astrometry_image_path)
        time.sleep(1)
        print('')

    #If the header DOES have a HISTORY keyword, skip it, it has already been processed. 
    else:
        print('Astrometric processing already complete for {}, skipping.'.format(filename.name))
        print('')

def source_pixels_to_world(short_name, source_detect_image_path, force_output_path=''): 
    """Gets world coordinates of tracked sources (target + references) in the source_detect_image. 

    :param short_name: short name for the target
    :type short_name: str
    :param source_detect_image_path: path to the source_detect_image
    :type source_detect_image_path: pathlib.PosixPath
    :param force_output_path: if you want to manually set an output directory, specify the top-level here (i.e. the directory containing Logs, Objects, etc.), defaults to ''
    :type force_output_path: str, optional
    :return: updates the csv of target/reference source_detect_centroids with world coordinates
    :rtype: csv
    """

    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pat.utils.pines_dir_check()

    sources_csv_path = pines_path/('Objects/'+short_name+'/sources/target_and_references_source_detection.csv')
    sources_df = pd.read_csv(sources_csv_path)

    #Get the WCS information. 
    source_detect_image = fits.open(source_detect_image_path)
    w = WCS(source_detect_image[0].header)

    #Get the world coordinates of all sources
    source_ras = np.zeros(len(sources_df))
    source_decs = np.zeros(len(sources_df))
    for i in range(len(sources_df)):
        source_x = sources_df['Source Detect X'][i]
        source_y = sources_df['Source Detect Y'][i]
        sky = w.pixel_to_world(source_x, source_y)
        source_ras[i] = sky.ra.value
        source_decs[i] = sky.dec.value
    
    #Add world coordinates to source detection df. 
    sources_df['Source Detect RA'] = source_ras
    sources_df['Source Detect Dec'] = source_decs

    #Write out to update the source detection csv. 
    sources_df.to_csv(sources_csv_path, index=0)
    return sources_df
   
def upload_file(apikey, filename, header):
    """Uploads files to astrometry.net. 

    :param apikey: api key of your account on astrometry.net
    :type apikey: str
    :param filename: path to the file to upload
    :type filename: pathlib.PosixPath
    :param header: header of the file you want to upload
    :type header: astropy.io fits header 
    """
    astrometry = Client()

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
            
