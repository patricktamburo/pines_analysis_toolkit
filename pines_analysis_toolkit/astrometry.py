from __future__ import print_function

import pines_analysis_toolkit as pat

from astropy.convolution import interpolate_replace_nans, Gaussian2DKernel
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u

import matplotlib.pyplot as plt

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
from email.mime.application import MIMEApplication
from email.encoders import encode_noop

import json

from astroquery.vizier import Vizier
from astroquery.gaia import Gaia

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
                 apiurl=default_url):
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
            args.update({'session': self.session})
        #print('Python:', args)
        json = python2json(args)
        #print('Sending json:', json)
        url = self.get_url(service)
        #print('Sending to URL:', url)

        # If we're sending a file, format a multipart/form-data
        if file_args is not None:
            import random
            boundary_key = ''.join(
                [random.choice('0123456789') for i in range(19)])
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
        args = {'apikey': apikey}
        result = self.send_request('login', args)
        sess = result.get('session')
        #print('Got session:', sess)
        if not sess:
            raise RequestError('no session in result')
        self.session = sess

    def _get_upload_args(self, **kwargs):
        args = {}
        for key, default, typ in [('allow_commercial_use', 'd', str),
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
                                  ('parity', None, int),
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
        result = self.send_request('submission_images', {'subid': subid})
        return result.get('image_ids')

    def overlay_plot(self, service, outfn, wcsfn, wcsext=0):
        from astrometry.util import util as anutil
        wcs = anutil.Tan(wcsfn, wcsext)
        params = dict(crval1=wcs.crval[0], crval2=wcs.crval[1],
                      crpix1=wcs.crpix[0], crpix2=wcs.crpix[1],
                      cd11=wcs.cd[0], cd12=wcs.cd[1],
                      cd21=wcs.cd[2], cd22=wcs.cd[3],
                      imagew=wcs.imagew, imageh=wcs.imageh)
        result = self.send_request(service, {'wcs': params})
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

    def annotate_data(self, job_id):
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

def gaia_cmd(short_name, sources, catalog='eDR3', plot=True, force_output_path=''):
    """Queries Gaia for sources and creates a CMD for them.
        Adds Gaia M_G and BP-RP to the source csv file


    :param short_name: the target's short name
    :type short_name: str
    :param sources: dataframe of sources, output from ref_star_chooser
    :type sources: pandas DataFrame
    :param catalog: which Gaia data release to query, either 'DR2' or 'eDR3', defaults to 'eDR3'
    :type catalog: str, optional
    :param plot: whether or not to save a plot of the CMD to the sources directory, defaults to True
    :type plot: bool, optional
    :param force_output_path: user-chosen top-level path for analysis if you don't want to use the default ~/Documents/PINES_analysis_toolkit/ folder, defaults to ''
    :type force_output_path: str, optional
    :raises ValueError: if catalog != 'DR2' | 'eDR3'
    :return: dataframe with new Gaia information added. Also updates the source csv.
    :rtype: pandas DataFrame
    """

    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pat.utils.pines_dir_check()

    # Set the Gaia data release to query.
    if catalog == 'DR2':
        Gaia.MAIN_GAIA_TABLE = "gaiadr2.gaia_source"
    elif catalog == 'eDR3':
        Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_source"
    else:
        raise ValueError('catalog must be DR2 or eDR3.')

    num_s = len(sources)
    bp_rp = np.zeros(num_s)  # Gaia Bp - Rp color
    p_mas = np.zeros(num_s)  # Parallax in mas
    M_G = np.zeros(num_s)  # Absolute Gaia Rp

    print('Querying sources in Gaia {}.'.format(catalog))

    for i in range(num_s):
        ra = sources['Source Detect RA'][i]
        dec = sources['Source Detect Dec'][i]
        c = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))

        # Query Gaia
        #query_gaia = Gaia.cone_search(c, radius=5*u.arcsec)
        #result_gaia = query_gaia.get_results()

        result_gaia = Gaia.query_object_async(coordinate=c, width=u.Quantity(
            5, u.arcsec), height=u.Quantity(5, u.arcsec))

        # If no sources at this position, return NaNs.
        if len(result_gaia) == 0:
            bp_rp[i] = np.nan
            p_mas[i] = np.nan
            M_G[i] = np.nan
            print('No source found in Gaia {} at ({:1.4f},{:1.4f}). Returning NaNs.'.format(
                catalog, ra, dec))
        # If one source found, return values
        else:
            bp_rp[i] = result_gaia['bp_rp'][0]
            m_G = result_gaia['phot_g_mean_mag'][0]
            p_mas[i] = result_gaia['parallax'][0]
            dist_pc = 1/(p_mas[i]/1000)
            M_G[i] = m_G + 5 - 5*np.log10(dist_pc)

    if plot:
        plt.figure(figsize=(6, 9))
        plt.scatter(bp_rp, M_G)
        plt.title(sources['Name'][0]+' field CMD', fontsize=18)
        plt.xlabel('G$_{BP}$ - G$_{RP}$', fontsize=16)
        plt.ylabel('M$_G$', fontsize=16)
        plt.xlim(-1, 5)
        plt.ylim(-5, 16)
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(pines_path/('Objects/'+short_name +
                    '/sources/gaia_cmd.png'), dpi=300)

    sources['p (mas)'] = p_mas
    sources['M_G'] = M_G
    sources['bp_rp'] = bp_rp

    sources = source_spts(sources)

    sources.to_csv(pines_path/('Objects/'+short_name +
                   '/sources/target_and_references_source_detection.csv'), index=0, na_rep='NaN')

    return sources

def source_spts(sources):
    """Uses a table from Eric Mamajek to estimate reference SpTs from their Gaia Bp-Rp colors.
    Uses a Teff-SpT relation from Faherty et al. (2016) to generate temperature estimates for the target.

    :param sources: dataframe of sources with Gaia Bp-Rp column
    :type sources: pandas DataFrame
    :return: dataframe with source SpTs added.
    :rtype: pandas DataFrame
    """

    pines_path = pat.utils.pines_dir_check()
    mamajek_path = pines_path/('Misc/mamajek_spts.txt')
    mamajek_df = pd.read_csv(mamajek_path, sep=r"[ ]{1,}")
    mamajek_spectral_types = np.array(mamajek_df['#SpT'])
    mamajek_teffs = np.array(mamajek_df['Teff'])

    mamajek_bp_rp = np.array(mamajek_df['Bp-Rp'])
    mamajek_bp_rp[mamajek_bp_rp == '...'] = np.nan
    mamajek_bp_rp = np.array([float(i) for i in mamajek_bp_rp])
   
   
    #Set up some lines, approximately parallel to the main sequence, that sources need to lie between in M_G vs. BP-RP to be considered dwarfs. Otherwise, discard them for being WDs or Giants.
    
    bp_rp_x = np.linspace(-0.5, 5.5, 1000)
    slope = 3
    lower_line = slope*bp_rp_x + 5 #Approximate WD separator
    upper_line = slope*bp_rp_x - 3 #Approximate giant separator
    
    source_bp_rps = np.array(sources['bp_rp'])
    source_M_Gs = np.array(sources['M_G'])
    source_spts = []
    source_teffs = []
    for i in range(len(source_bp_rps)):
        source_bp_rp = source_bp_rps[i]
        source_M_G = source_M_Gs[i]
        if not np.isnan(source_bp_rp):
            wd_giant_check_ind = np.where(abs(bp_rp_x-source_bp_rp) == np.min(abs(bp_rp_x-source_bp_rp)))[0][0]
        if (not np.isnan(source_bp_rp)) and (source_bp_rp > -0.120) and (source_bp_rp < 4.8) and (source_M_G < lower_line[wd_giant_check_ind]) and (source_M_G > upper_line[wd_giant_check_ind]):
            closest_spt_ind = np.where(abs(
                mamajek_bp_rp-source_bp_rp) == np.nanmin(abs(mamajek_bp_rp-source_bp_rp)))[0][0]
            source_spts.append(mamajek_spectral_types[closest_spt_ind])
            source_teffs.append(mamajek_teffs[closest_spt_ind])
        else:
            if i != 0:
                #TODO: Have to update to actually toss references identified as WDs or giants.
                breakpoint()
            source_spts.append(np.nan)
            source_teffs.append(np.nan)

    sources['SpT'] = source_spts
    sources['Teff'] = source_teffs
    
    #Add info for the target, which probably wasn't found in Gaia
    pines_sample_path = pines_path/('Misc/PINES Sample.xlsx')
    pines_sample_df = pd.read_excel(pines_sample_path)
    sample_entry = np.where(np.array(pines_sample_df['Short Name'] == 'J'+sources['Name'][0].split(' ')[1]))[0][0]
    target_spt = pines_sample_df['SpT'][sample_entry]
    sources['SpT'][0] = target_spt
    target_temp = faherty_2016_spt_teff_relation(target_spt)
    sources['Teff'][0] = np.round(target_temp)
    return sources

def faherty_2016_spt_teff_relation(spt):
    """Calculates a temperature given a spectral type using a relation from Faherty et al. (2016)
    https://iopscience.iop.org/article/10.3847/0067-0049/225/1/10/pdf#%FE%FF%00b%00m%00_%00a%00p%00j%00s%00a%00a%002%005%009%00d%00t%001%009


    :param spt: string specifying the spectral type (e.g., 'L5')
    :type spt: str
    :return: temperature of the spectral type
    :rtype: float
    """
    spectral_class = spt[0]
    if spectral_class == 'M':
        num = float(spt[1:])
    elif spectral_class == 'L':
        num = float(spt[1:]) + 10
    elif spectral_class == 'T':
        num = float(spt[1:]) + 20
    
    #Coefficients from T_eff FLD relation, Table 19 of Faherty + 2016. 
    c0 = 4.747e3
    c1 = -7.005e2
    c2 = 1.155e2
    c3 = -1.191e1
    c4 = 6.318e-1
    c5 = -1.606e-2
    c6 = 1.546e-4

    temp = c0*num**0 + c1*num**1 + c2*num**2 + c3*num**3 + c4*num**4 + c5*num**5 + c6*num**6

    return temp

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
        print(short_name+', '+filename.name+', image ' +
              str(i+1)+' of '+str(len(reduced_files)))

        # If the header does not have a HISTORY keyword (which is added by astrometry.net), process it.
        header = fits.open(filename)[0].header

        if 'HISTORY' not in header:
            # Read in the image data, interpolate NaNs, and save to a temporary fits file.
            # Astrometry.net does not work with NaNs in images.
            original_image = fits.open(filename)[0].data
            image = interpolate_replace_nans(original_image, kernel=kernel)
            temp_filename = filename.parent / \
                (filename.name.split('.fits')[0]+'_temp.fits')
            hdu = fits.PrimaryHDU(image, header=header)
            hdu.writeto(temp_filename, overwrite=True)

            # Upload the image with interpolated NaNs to astrometry.net. It will solve and download to a new image.
            pat.astrometry.upload_to_astrometry.upload_file(
                api_key, temp_filename, header)

            # Try to donwload the solved image and open it.
            try:
                # Grab the header of the astrometry.net solution image, and the original image data.
                astrometry_image_path = filename.parent / \
                    (temp_filename.name.split('.fits')[0]+'_new_image.fits')
                wcs_header = fits.open(astrometry_image_path)[0].header
                wcs_hdu = fits.PrimaryHDU(original_image, header=wcs_header)

                # Save the original image data with the new wcs header.
                output_filename = filename
                wcs_hdu.writeto(output_filename, overwrite=True)

            # If the try clause didn't work, that's because the processing on astrometry.net failed.
            except:
                print('Astrometry solution failed!')
                # Add HISTORY keyword to header so it doesn't get processed in the future.
                header['HISTORY'] = 'Astrometric solution failed.'
                fits.writeto(filename, original_image, header, overwrite=True)

            # Delete temporary files.
            os.remove(temp_filename)
            os.remove(astrometry_image_path)
            time.sleep(1)
            t2 = time.time()
            times.append(t2-t1)
            print('Average completion time: {:3.1f} s.'.format(np.mean(times)))
            print('')

        # If the header DOES have a HISTORY keyword, skip it, it has already been processed.
        else:
            print('Astrometric processing already complete for {}, skipping.'.format(
                filename.name))
            print('')


def source_detect_astrometry(short_name, api_key, filename, filter, force_output_path=''):
    """ Uploads a single reduced image for a target to astrometry.net, downloads solution image, and updates the image header with the astrometry.net wcs. 
    
    :param short_name: the short name of the target
    :type short_name: str
    :param api_key: the api key of your account on astrometry.net
    :type api_key: str
    :param filename: name of the image you want to process
    :type filename: str
    :param force_output_path: if you want to manually set an output directory, specify the top-level here (i.e. the directory containing Logs, Objects, etc.), defaults to ''
    :type force_output_path: str, optional
    :raises RuntimeError: if astrometry solution fails
    """
    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pat.utils.pines_dir_check()
    
    image_path = pines_path/('Objects/'+short_name+'/reduced/'+filter+'/'+filename)
    
    kernel = Gaussian2DKernel(x_stddev=0.5)

    # If the header does not have a HISTORY keyword (which is added by astrometry.net), process it.
    header = fits.open(image_path)[0].header

    if 'HISTORY' not in header:
        print('Uploading {} to astrometry.net.'.format(image_path.name))
        # Read in the image data, interpolate NaNs, and save to a temporary fits file.
        # Astrometry.net does not work with NaNs in images.
        original_image = fits.open(image_path)[0].data
        image = interpolate_replace_nans(original_image, kernel=kernel)
        temp_filename = image_path.parent / \
            (image_path.name.split('.fits')[0]+'_temp.fits')
        hdu = fits.PrimaryHDU(image, header=header)
        hdu.writeto(temp_filename, overwrite=True)

        # Upload the image with interpolated NaNs to astrometry.net. It will solve and download to a new image.
        pat.astrometry.upload_file(api_key, temp_filename, header)

        # Try to donwload the solved image and open it.
        try:
            # Grab the header of the astrometry.net solution image, and the original image data.
            astrometry_image_path = image_path.parent / \
                (temp_filename.name.split('.fits')[0]+'_new_image.fits')
            wcs_header = fits.open(astrometry_image_path)[0].header
            wcs_hdu = fits.PrimaryHDU(original_image, header=wcs_header)

            # Save the original image data with the new wcs header.
            output_filename = image_path
            wcs_hdu.writeto(output_filename, overwrite=True)

        # If the try clause didn't work, that's because the processing on astrometry.net failed.
        except:
            raise RuntimeError(
                'Astrometry solution failed! Use a different source detect image for choosing reference stars.')

        # Delete temporary files.
        os.remove(temp_filename)
        os.remove(astrometry_image_path)
        time.sleep(1)
        print('')

    # If the header DOES have a HISTORY keyword, skip it, it has already been processed.
    else:
        print('Astrometric processing already complete for {}, skipping.'.format(
            image_path.name))
        print('')


def source_pixels_to_world(short_name, wcs_filename, filter, force_output_path=''):
    """Gets world coordinates of tracked sources (target + references) in the source_detect_image. 

    :param short_name: short name for the target
    :type short_name: str
    :param wcs_filename: name of the file with wcs information from source_detect_astrometry processing
    :type wcs_filename: str
    :param force_output_path: if you want to manually set an output directory, specify the top-level here (i.e. the directory containing Logs, Objects, etc.), defaults to ''
    :type force_output_path: str, optional
    :return: updates the csv of target/reference source_detect_centroids with world coordinates
    :rtype: csv
    """

    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pat.utils.pines_dir_check()

    sources_csv_path = pines_path / \
        ('Objects/'+short_name+'/sources/target_and_references_source_detection.csv')
    sources_df = pd.read_csv(sources_csv_path)

    # Get the WCS information.
    source_detect_image_path = pines_path/('Objects/'+short_name+'/reduced/'+filter+'/'+wcs_filename)
    source_detect_image = fits.open(source_detect_image_path)
    w = WCS(source_detect_image[0].header)

    # Get the world coordinates of all sources
    source_ras = np.zeros(len(sources_df))
    source_decs = np.zeros(len(sources_df))
    for i in range(len(sources_df)):
        source_x = sources_df['Source Detect X'][i]
        source_y = sources_df['Source Detect Y'][i]
        sky = w.pixel_to_world(source_x, source_y)
        source_ras[i] = sky.ra.value
        source_decs[i] = sky.dec.value

    # Add world coordinates to source detection df.
    sources_df['Source Detect RA'] = source_ras
    sources_df['Source Detect Dec'] = source_decs

    # Write out to update the source detection csv.
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
    tel_ra_deg = int(tel_ra_str[0])*15 + int(tel_ra_str[1]) * \
        (15/60) + float(tel_ra_str[2])*(15/3600)

    tel_dec_str = header['TELDEC'].split(':')

    if tel_dec_str[0][0] == '+':
        tel_dec_deg = int(
            tel_dec_str[0]) + int(tel_dec_str[1])/60 + float(tel_dec_str[2])/3600
    elif tel_dec_str[0][0] == '-':
        tel_dec_deg = int(
            tel_dec_str[0]) - int(tel_dec_str[1])/60 - float(tel_dec_str[2])/3600

    sub_id = astrometry.upload(filename, scale_units='arcsecperpix', scale_lower=0.575,
                               scale_upper=0.585, center_ra=tel_ra_deg, center_dec=tel_dec_deg, radius=0.116)['subid']
    print('submission_id: ', str(sub_id))

    while astrometry.sub_status(sub_id)['jobs'] == [] or astrometry.sub_status(sub_id)['jobs'] == [None]:
        time.sleep(1)

    job_id = astrometry.sub_status(sub_id)['jobs'][0]
    print('job_id: ', astrometry.sub_status(sub_id)['jobs'])

    while astrometry.job_status(job_id) == 'solving':
        time.sleep(1)

    outfile = filename.parent / \
        (filename.name.split('.fits')[0]+'_new_image.fits')

    url = 'http://nova.astrometry.net/new_fits_file/' + str(job_id)

    breakpoint()
    r = requests.get(url, allow_redirects=True)
    open(outfile, 'wb').write(r.content)
