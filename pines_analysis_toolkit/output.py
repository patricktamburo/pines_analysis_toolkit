from pines_analysis_toolkit.utils import pines_dir_check, short_name_creator, pines_login

from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.pagesizes import LETTER
from reportlab.lib.utils import ImageReader

import pdb
import numpy as np
from natsort import natsorted
from datetime import date, datetime
import os
import shutil
from numpy.lib.utils import source
from pathlib import Path


def dv_report(target, pat_version='1.0', force_output_path=''):

    def page_header():
        today = date.today()
        date_str = today.strftime('%b %d, %Y')
        time = datetime.now()
        time_str = time.strftime('%H:%M')
        canvas.setFont('Times-Italic', 8)  # Font for captions
        canvas.setFillColor('Grey')
        header_text = short_name+' DV Report, Compiled on '+date_str+' at ' + \
            time_str+' with PINES Analysis Toolkit V. '+str(pat_version)
        canvas.drawCentredString(612*0.5, 792-15, header_text)

    def page_footer():
        canvas.setFont('Times-Italic', 8)  # Font for captions
        canvas.setFillColor('Grey')
        footer_text = str(page_num)
        canvas.drawCentredString(612*0.95, 792*0.02, footer_text)

    def title_page():
        page_header()
        global fig_num, page_num
        canvas.setFont('Times-Roman', 20)  # Font for captions
        canvas.setFillColor('Black')
        title_text = 'PINES DV Report'
        canvas.drawCentredString(fig_x*0.5, fig_y*0.9, title_text)
        object_text = target+' ('+short_name+')'
        canvas.setFont('Times-Roman', 16)  # Font for captions
        canvas.drawCentredString(fig_x*0.5, fig_y*0.87, object_text)
        # Add the source image.
        source_image_path = source_path/('target_and_refs.png')
        img = ImageReader(source_image_path)
        y = 50
        w = fig_x  # All the way across the page
        h = w * source_aspect_ratio
        canvas.drawImage(img, x, y, x+w, h)
        canvas.setFont('Times-Roman', 12)  # Font for captions
        caption_text = 'Figure {}: Target and references.'.format(fig_num)
        canvas.drawCentredString(x+w*0.5, 30, caption_text)
        fig_num += 1
        page_footer()

        canvas.showPage()
        page_num += 1

    def night_page():
        global fig_num, page_num
        # Get all nightly-normalized target lightcurve plots in the analysis directory.
        nightly_glob = natsorted(
            np.array([x for x in analysis_path.glob('*nightly*target_flux.png')]))

        if os.path.exists(analysis_path/('optimal_aperture.txt')):
            with open(analysis_path/('optimal_aperture.txt'), 'r') as f:
                lines = f.readlines()
                best_ap = lines[0].split('  ')[1].split('\n')[0]
                ap_type = best_ap.split('_')[1]
        best_aper_phot_path = analysis_path/('aper_phot_analysis/'+best_ap)

        nightly_glob = np.array(
            natsorted([x for x in best_aper_phot_path.glob('*night*.png')]))
        # Sort plots in order we want them in the report: arget lc, raw flux, normalized flux
        nightly_glob = np.array(
            [nightly_glob[2], nightly_glob[1], nightly_glob[0]])

        if ap_type == 'fixed':
            rad = best_ap.split('_')[0]
            cap_ender = 'radius = {} pixels'.format(rad)
        elif ap_type == 'variable':
            fact = best_ap.split('_')[0]
            cap_ender = 'multiplicative seeing factor = {}'.format(fact)
            page_footer

        caption_texts = ['Best nightly-normalized corrected target flux, {} aperture photometry, {}.'.format(ap_type, cap_ender),
                         'Raw flux, {} aperture photometry, {}.'.format(
                             ap_type, cap_ender),
                         'Nightly-normalized flux, {} aperture photometry, {}.'.format(ap_type, cap_ender)]

        for i in range(len(nightly_glob)):
            if i == 0:
                page_header()
            canvas.setFont('Times-Roman', 12)  # Font for captions
            canvas.setFillColor('Black')
            image_path = str(nightly_glob[i])

            img = ImageReader(image_path)
            y = (3-(i+1))*fig_y/3 + 50
            w = fig_x  # All the way across the page
            h = w * lc_aspect_ratio
            canvas.drawImage(img, x, y, x+w, h)
            caption_text = 'Figure {}: '.format(fig_num)+caption_texts[i]
            canvas.drawCentredString(x+w*0.5, y-h*0.1, caption_text)
            fig_num += 1

            if i == len(nightly_glob) - 1:
                page_footer()

        canvas.showPage()
        page_num += 1

    # def global_page():
    #     global fig_num, page_num

    #     if os.path.exists(analysis_path/('optimal_aperture.txt')):
    #         with open(analysis_path/('optimal_aperture.txt'), 'r') as f:
    #             lines = f.readlines()
    #             best_ap = lines[1].split(' ')[1]
    #             ap_type = best_ap.split('_')[1]
    #     best_aper_phot_path = analysis_path/('aper_phot_analysis/'+best_ap)

    #     global_glob = np.array(natsorted([x for x in best_aper_phot_path.glob('*global*.png')]))
    #     #Sort plots in order we want them in the report: arget lc, raw flux, normalized flux
    #     global_glob = np.array([global_glob[2], global_glob[1], global_glob[0]])

    #     if ap_type == 'fixed':
    #         rad = best_ap.split('_')[0]
    #         cap_ender = 'radius = {} pixels'.format(rad)
    #     elif ap_type == 'variable':
    #         fact = best_ap.split('_')[0]
    #         cap_ender = 'multiplicative seeing factor = {}'.format(fact)
    #         page_footer

    #     caption_texts = ['Best globally-normalized corrected target flux, {} aperture photometry, {}.'.format(ap_type, cap_ender),
    #                     'Raw flux, {} aperture photometry, {}.'.format(ap_type, cap_ender),
    #                     'Globally-normalized flux, {} aperture photometry, {}.'.format(ap_type, cap_ender)]

    #     for i in range(len(global_glob)):
    #         if i == 0:
    #             page_header()
    #         canvas.setFont('Times-Roman', 12) #Font for captions
    #         canvas.setFillColor('Black')
    #         image_path = str(global_glob[i])
    #         aperture_radius = image_path.split('=')[1].split('_')[0]

    #         img = ImageReader(image_path)
    #         y = (3-(i+1))*fig_y/3 + 50
    #         w = fig_x #All the way across the page
    #         h = w * lc_aspect_ratio
    #         canvas.drawImage(img, x, y, x+w, h)
    #         caption_text = 'Figure {}: '.format(fig_num)+caption_texts[i]
    #         canvas.drawCentredString(x+w*0.5,y-h*0.1, caption_text)
    #         fig_num += 1

    #         if i == len(global_glob) - 1:
    #             page_footer()

    #     canvas.showPage()
    #     page_num += 1

    def centroid_page():
        global fig_num, page_num
        page_header()

        canvas.setFont('Times-Roman', 12)  # Font for captions
        canvas.setFillColor('Black')
        cutout_position_path = diagnostic_plot_path / \
            (short_name+'_cutout_positions.png')
        img = ImageReader(cutout_position_path)
        y = 2*fig_y/3 - 100
        w = fig_x  # All the way across the page
        h = w * centroid_aspect_ratio
        canvas.drawImage(img, x, y, x+w, h)
        caption_text = 'Figure {}: '.format(
            fig_num)+'Cutout image positions for all sources.'
        canvas.drawCentredString(x+w*0.5, y-h*0.05, caption_text)
        fig_num += 1

        image_position_path = diagnostic_plot_path / \
            (short_name+'_image_positions.png')
        img = ImageReader(image_position_path)
        y = 1*fig_y/3 - 200
        w = fig_x  # All the way across the page
        h = w * centroid_aspect_ratio
        canvas.drawImage(img, x, y, x+w, h)
        caption_text = 'Figure {}: '.format(
            fig_num)+'Absolute image positions for the target.'
        canvas.drawCentredString(x+w*0.5, y-h*0.05, caption_text)
        fig_num += 1

        page_footer()
        canvas.showPage()
        page_num += 1

    def diagnostics_page():
        global fig_num, page_num
        page_header()

        canvas.setFont('Times-Roman', 12)  # Font for captions
        canvas.setFillColor('Black')
        seeing_path = diagnostic_plot_path/(short_name+'_seeing.png')
        img = ImageReader(seeing_path)
        y = 2*fig_y/3 + 50
        w = fig_x  # All the way across the page
        h = w * lc_aspect_ratio
        canvas.drawImage(img, x, y, x+w, h)
        caption_text = 'Figure {}: '.format(
            fig_num)+'Seeing measurements (in arcsec).'
        canvas.drawCentredString(x+w*0.5, y-h*0.1, caption_text)
        fig_num += 1

        background_path = diagnostic_plot_path/(short_name+'_backgrounds.png')
        img = ImageReader(background_path)
        y = 1*fig_y/3 + 50
        w = fig_x  # All the way across the page
        h = w * lc_aspect_ratio
        canvas.drawImage(img, x, y, x+w, h)
        caption_text = 'Figure {}: '.format(
            fig_num)+'Background measurements (in ADU). Non-linear effects begin near 4000 ADU.'
        canvas.drawCentredString(x+w*0.5, y-h*0.1, caption_text)
        fig_num += 1

        airmass_path = diagnostic_plot_path/(short_name+'_airmasses.png')
        img = ImageReader(airmass_path)
        y = 0*fig_y/3 + 50
        w = fig_x  # All the way across the page
        h = w * lc_aspect_ratio
        canvas.drawImage(img, x, y, x+w, h)
        caption_text = 'Figure {}: '.format(fig_num)+'Airmass measurements.'
        canvas.drawCentredString(x+w*0.5, y-h*0.1, caption_text)
        fig_num += 1

        page_footer()

        canvas.showPage()

    def corr_refs_pages():
        global fig_num, page_num

        if os.path.exists(analysis_path/('optimal_aperture.txt')):
            with open(analysis_path/('optimal_aperture.txt'), 'r') as f:
                lines = f.readlines()
                best_ap = lines[0].split('  ')[1].split('\n')[0]
                ap_type = best_ap.split('_')[1]
        best_aper_phot_path = analysis_path / \
            ('aper_phot_analysis/'+best_ap+'/corr_ref_plots/')
        corr_glob = np.array(
            natsorted(list([x for x in best_aper_phot_path.glob('*.png')])))

        count = 0
        for i in range(len(corr_glob)):
            if count % 3 == 0:
                page_header()
            canvas.setFont('Times-Roman', 12)  # Font for captions
            canvas.setFillColor('Black')
            image_path = str(corr_glob[i])

            img = ImageReader(image_path)
            y = (3-((count % 3)+1))*fig_y/3 + 50
            w = fig_x  # All the way across the page
            h = w * lc_aspect_ratio
            canvas.drawImage(img, x, y, x+w, h)
            if i == 0:
                caption_text = 'Figure {}: '.format(
                    fig_num)+image_path.split('/')[-1].split('_')[0]+' corrected flux.'
            else:
                caption_text = 'Figure {}: '.format(
                    fig_num)+'Reference '+str(i)+' corrected flux.'

            canvas.drawCentredString(x+w*0.5, y-h*0.1, caption_text)
            fig_num += 1

            count += 1
            if count % 3 == 0:
                page_footer()
                canvas.showPage()
                page_num += 1

    # Set up pathing.
    if force_output_path != '':
        pines_path = force_output_path
    else:
        pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    print('Generating PINES DV report for {}.'.format(short_name))
    output_filename = (short_name+'_dv_report.pdf').replace(' ', '')
    target_path = pines_path/('Objects/'+short_name)
    source_path = target_path/('sources/')
    analysis_path = target_path/('analysis/')
    diagnostic_plot_path = analysis_path/('diagnostic_plots/')
    output_path = target_path/('output/'+output_filename)

    # Set up the DV report
    canvas = Canvas(str(output_path), pagesize=LETTER)
    global fig_num, page_num, fig_x, fig_y, source_aspect_ratio, lc_aspect_ratio, centroid_aspect_ratio, x
    fig_num = 1  # Initialize
    page_num = 1
    fig_x = 612
    fig_y = 792  # Dimensions (in pixels?) for 8.5x11 page
    source_aspect_ratio = 9/10
    lc_aspect_ratio = 5/17
    centroid_aspect_ratio = 27/51
    x = 0  # Start all figures at the left margin

    # ------------------------------------------------------------------------------------------------------------------------
    # PAGE 1: Title and Source Detection Image.
    # Make the title page.
    title_page()

    # ------------------------------------------------------------------------------------------------------------------------
    # PAGE 2: Best nightly-normalized lightcurve, raw flux, and normalized flux.
    night_page()

    # ------------------------------------------------------------------------------------------------------------------------
    # #PAGE 3: Best globally-normalized lightcurve, raw flux, and normalized flux.
    # global_page()

    # ------------------------------------------------------------------------------------------------------------------------

    # PAGE 3: Centroid diagnostic plots.
    centroid_page()

    # ------------------------------------------------------------------------------------------------------------------------
    # PAGE 4: Diagnostic plots.
    diagnostics_page()

    # ------------------------------------------------------------------------------------------------------------------------
    # PAGES 5-N: Plots of corrected target flux.
    corr_refs_pages()

    canvas.save()

    return


def output_wrangler(target):
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    object_path = pines_path/('Objects/'+short_name)
    sources_path = object_path/('sources')
    aper_phot_path = object_path/('aper_phot')
    analysis_path = object_path/('analysis')
    output_path = object_path/('output')

    # Copy souce detect image.
    src = sources_path/('target_and_refs.png')
    dest = output_path/('target_and_refs.png')
    shutil.copyfile(src, dest)

    # Copy image centroids.
    src = sources_path/('target_and_references_centroids.csv')
    dest = output_path/('target_and_references_centroids.csv')
    shutil.copyfile(src, dest)

    # Copy best aperture raw photometry.
    best_ap_path = analysis_path/('optimal_aperture.txt')
    with open(best_ap_path, 'r') as f:
        best_ap = f.readlines()[0]
    src = aper_phot_path/(short_name+'_aper_phot_'+best_ap+'_pix_radius.csv')
    dest = output_path/(short_name+'_aper_phot_'+best_ap+'_pix_radius.csv')
    shutil.copyfile(src, dest)

    # Copy best aperture weighted lightcurve photometry.
    src = analysis_path/('aper_phot_analysis/'+best_ap+'/' +
                         short_name+'_weighted_lc_aper_phot_'+best_ap+'_pix.csv')
    dest = output_path / \
        (short_name+'_weighted_lc_aper_phot_'+best_ap+'_pix.csv')
    shutil.copyfile(src, dest)

    pdb.set_trace()
    return


def upload_output(sftp, target):
    ''' AUTHORS
            Patrick Tamburo, BU, August 2021
        PURPOSE
            Uploads output to the /data/lightcurves folder on the PINES server.
        INPUTS
            sftp (pysftp connection): sftp connection to the PINES server.
            target (str): Long name for the target. 
        OUTPUTS
            None.
        TODO:
            None.
    '''
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)

    # Get the path of the object's DV report on the local machine.
    dv_report_path = pines_path / \
        ('Objects/'+short_name+'/output/' +
         short_name.replace(' ', '')+'_dv_report.pdf')

    # Move the sftp connection to the object's lightcurves/ directory on the PINES server.
    sftp.chdir('/data/lightcurves/')
    # If the folder doesn't exist, make it.
    if not short_name in sftp.listdir():
        sftp.mkdir(short_name)
    sftp.chdir(short_name)

    # Upload the DV report
    if os.path.exists(dv_report_path):
        print('Uploading {} DV report.'.format(short_name))
        upload_path = sftp.pwd+'/'+short_name.replace(' ', '')+'_dv_report.pdf'
        sftp.put(dv_report_path, upload_path)
    else:
        print('WARNING: No DV report found for {}.'.format(short_name))

    if not 'sources' in sftp.listdir():
        sftp.mkdir('sources')

    sftp.chdir('sources')

    # Upload source detection CSV.
    source_detection_path = pines_path / \
        ('Objects/'+short_name+'/sources/target_and_references_source_detection.csv')
    if os.path.exists(source_detection_path):
        print('Uploading {} source detection CSV.'.format(short_name))
        upload_path = sftp.pwd+'/target_and_references_source_detection.csv'
        sftp.put(source_detection_path, upload_path)
    else:
        print('WARNING: No source detection CSV found for {}.'.format(short_name))

    # Upload centroids CSV.
    centroid_path = pines_path / \
        ('Objects/'+short_name+'/sources/target_and_references_centroids.csv')
    if os.path.exists(centroid_path):
        print('Uploading {} centroid CSV.'.format(short_name))
        upload_path = sftp.pwd+'/target_and_references_centroids.csv'
        sftp.put(centroid_path, upload_path)
    else:
        print('WARNING: No centroid CSV found for {}.'.format(short_name))
    sftp.chdir('..')

    if not 'aper_phot' in sftp.listdir():
        sftp.mkdir('aper_phot')
    sftp.chdir('aper_phot')

    aper_phot_path = pines_path/('Objects/'+short_name+'/aper_phot/')
    aper_phot_files = np.array(
        natsorted([x for x in aper_phot_path.glob('*.csv')]))
    n_aper_phot_files = len(aper_phot_files)
    if n_aper_phot_files != 0:
        for i in range(n_aper_phot_files):
            if i == 0:
                print('Uploading {} aper_phot CSVs.'.format(short_name))
            aper_phot_path = aper_phot_files[i]
            upload_path = sftp.pwd+'/'+aper_phot_path.name
            sftp.put(aper_phot_path, upload_path)
    else:
        print('WARNING: No aper_phot CSVs found for {}.'.format(short_name))
    sftp.chdir('..')

    if not 'analysis' in sftp.listdir():
        sftp.mkdir('analysis')
    sftp.chdir('analysis')

    root = pines_path/('Objects/'+short_name+'/analysis')
    f = []
    for path, subdirs, files in os.walk(root):
        for name in files:
            if 'weighted_lc' in name:
                f.append(Path(os.path.join(path, name)))
    n_weighted_lc_files = len(f)
    if n_aper_phot_files != 0:
        for i in range(n_weighted_lc_files):
            if i == 0:
                print('Uploading {} weighted_lc CSVs.'.format(short_name))
            weighted_lc_path = f[i]
            upload_path = sftp.pwd+'/'+weighted_lc_path.name
            sftp.put(weighted_lc_path, upload_path)
    else:
        print('WARNING: No weighted_lc CSVs found for {}.'.format(short_name))

    print('')
    return
