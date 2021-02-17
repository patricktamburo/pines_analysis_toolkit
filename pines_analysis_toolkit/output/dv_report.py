from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.pagesizes import LETTER
from reportlab.lib.utils import ImageReader
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
import pdb 
import numpy as np 
from natsort import natsorted
from datetime import date
import os 

def dv_report(target, pat_version='1.0'):
    
    def page_header():
        today = date.today()
        date_str = today.strftime('%b %d, %Y')
        canvas.setFont('Times-Italic', 8) #Font for captions
        canvas.setFillColor('Grey')
        header_text = short_name+' DV Report, Compiled on '+date_str+' with PINES Analysis Toolkit V. '+str(pat_version)
        canvas.drawCentredString(612*0.5, 792-15, header_text)

    def page_footer():
        canvas.setFont('Times-Italic', 8) #Font for captions
        canvas.setFillColor('Grey')
        footer_text = str(page_num)
        canvas.drawCentredString(612*0.95, 792*0.02, footer_text)

    def title_page():
        canvas.setFont('Times-Roman', 20) #Font for captions
        canvas.setFillColor('Black')
        title_text = 'PINES DV Report'
        canvas.drawCentredString(fig_x*0.5, fig_y*0.9, title_text)
        object_text = target+' ('+short_name+')'
        canvas.setFont('Times-Roman', 16) #Font for captions
        canvas.drawCentredString(fig_x*0.5, fig_y*0.87, object_text)

    #Set up pathing.
    pines_path = pines_dir_check()
    short_name = short_name_creator(target)
    print('Generating PINES DV report for {}.'.format(short_name))
    output_filename = (short_name+'_dv_report.pdf').replace(' ','')
    target_path = pines_path/('Objects/'+short_name)
    source_path = target_path/('sources/')
    analysis_path = target_path/('analysis/')
    diagnostic_plot_path = analysis_path/('diagnostic_plots/')
    output_path = target_path/('output/'+output_filename)

    #Set up the DV report
    canvas = Canvas(str(output_path), pagesize=LETTER)
    fig_num = 1 #Initialize 
    page_num = 1 
    fig_x = 612
    fig_y = 792 #Dimensions (in pixels?) for 8.5x11 page
    source_aspect_ratio = 9/10
    lc_aspect_ratio = 5/17
    centroid_aspect_ratio = 27/51

    #------------------------------------------------------------------------------------------------------------------------
    #PAGE 1: Title and Source Detection Image.
    #Make the title page. 
    page_header()
    title_page()

    #Add the source image. 
    source_image_path = source_path/('target_and_refs.png')
    img = ImageReader(source_image_path)
    x = 0 #Start at the left margin
    y = 50
    w = fig_x #All the way across the page
    h = w * source_aspect_ratio
    canvas.drawImage(img, x, y, x+w, h)
    canvas.setFont('Times-Roman', 12) #Font for captions
    caption_text = 'Figure {}: Target and references.'.format(fig_num)
    canvas.drawCentredString(x+w*0.5, 30, caption_text)
    fig_num += 1
    page_footer()

    canvas.showPage()
    page_num += 1

    #------------------------------------------------------------------------------------------------------------------------
    #PAGE 2: Best nightly-normalized lightcurve, raw flux, and normalized flux. 

    #Get all nightly-normalized target lightcurve plots in the analysis directory.
    nightly_glob = natsorted(np.array([x for x in analysis_path.glob('*nightly*target_flux.png')]))

    if os.path.exists(analysis_path/('optimal_aperture.txt')):
        with open(analysis_path/('optimal_aperture.txt'), 'r') as f:
            best_ap = f.readlines()[0]

    best_aper_phot_path = analysis_path/('aper_phot_analysis/'+best_ap)

    nightly_glob = np.array(natsorted([x for x in best_aper_phot_path.glob('*nightly*.png')]))
    nightly_glob = np.array([nightly_glob[3], nightly_glob[2], nightly_glob[1]])
    caption_texts = ['Best nightly-normalized corrected target flux, aperture radius of {} pixels.'.format(best_ap),
                    'Raw flux, aperture radius of {} pixels.'.format(best_ap),
                    'Nightly-normalized flux, aperture radius of {} pixels.'.format(best_ap)]
    
    for i in range(len(nightly_glob)):
        page_header()
        canvas.setFont('Times-Roman', 12) #Font for captions
        canvas.setFillColor('Black')
        image_path = str(nightly_glob[i])
        aperture_radius = image_path.split('=')[1].split('_')[0]
        
        img = ImageReader(image_path)
        y = (3-(i+1))*fig_y/3 + 30 
        w = fig_x #All the way across the page
        h = w * lc_aspect_ratio
        canvas.drawImage(img, x, y, x+w, h)
        caption_text = 'Figure {}: '.format(fig_num)+caption_texts[i]
        canvas.drawCentredString(x+w*0.5,y-h*0.1, caption_text)
        fig_num += 1

        if i == len(nightly_glob) - 1:
            page_footer()

    canvas.showPage()
    page_num += 1
    #------------------------------------------------------------------------------------------------------------------------
    
    #PAGE 3: Centroid diagnostic plots. 
    page_header()

    canvas.setFont('Times-Roman', 12) #Font for captions
    canvas.setFillColor('Black')
    cutout_position_path = diagnostic_plot_path/(short_name+'_cutout_positions.png')
    img = ImageReader(cutout_position_path)
    y = 2*fig_y/3 - 100
    w = fig_x #All the way across the page
    h = w * centroid_aspect_ratio
    canvas.drawImage(img, x, y, x+w, h)
    caption_text = 'Figure {}: '.format(fig_num)+'Cutout image positions for all sources.'
    canvas.drawCentredString(x+w*0.5,y-h*0.05, caption_text)
    fig_num += 1

    image_position_path = diagnostic_plot_path/(short_name+'_image_positions.png')
    img = ImageReader(image_position_path)
    y = 1*fig_y/3 - 200
    w = fig_x #All the way across the page
    h = w * centroid_aspect_ratio
    canvas.drawImage(img, x, y, x+w, h)
    caption_text = 'Figure {}: '.format(fig_num)+'Absolute image positions for the target.'
    canvas.drawCentredString(x+w*0.5,y-h*0.05, caption_text)
    fig_num += 1

    page_footer()
    canvas.showPage()
    page_num += 1

    #------------------------------------------------------------------------------------------------------------------------
    
    #PAGE 4: Ancillary diagnostic plots.
    page_header()

    canvas.setFont('Times-Roman', 12) #Font for captions
    canvas.setFillColor('Black')
    seeing_path = diagnostic_plot_path/(short_name+'_seeing.png')
    img = ImageReader(seeing_path)
    y = 2*fig_y/3 +30
    w = fig_x #All the way across the page
    h = w * lc_aspect_ratio
    canvas.drawImage(img, x, y, x+w, h)
    caption_text = 'Figure {}: '.format(fig_num)+'Seeing measurements (in arcsec).'
    canvas.drawCentredString(x+w*0.5,y-h*0.1, caption_text)
    fig_num += 1


    background_path = diagnostic_plot_path/(short_name+'_backgrounds.png')
    img = ImageReader(background_path)
    y = 1*fig_y/3 + 30
    w = fig_x #All the way across the page
    h = w * lc_aspect_ratio
    canvas.drawImage(img, x, y, x+w, h)
    caption_text = 'Figure {}: '.format(fig_num)+'Background measurements (in ADU).'
    canvas.drawCentredString(x+w*0.5,y-h*0.1, caption_text)
    fig_num += 1

    page_footer()


    canvas.showPage()

    canvas.save()
    
    pdb.set_trace()
    return

if __name__ == '__main__':
    target = '2MASS J00144919-0838207' 
    dv_report(target)