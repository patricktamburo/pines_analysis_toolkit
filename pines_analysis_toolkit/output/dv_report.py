from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.pagesizes import LETTER
from reportlab.lib.utils import ImageReader
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
import pdb 
import numpy as np 
from natsort import natsorted
from datetime import date, datetime
import os 

def dv_report(target, pat_version='1.0'):
    
    def page_header():
        today = date.today()
        date_str = today.strftime('%b %d, %Y')
        time = datetime.now()
        time_str = time.strftime('%H:%M')
        canvas.setFont('Times-Italic', 8) #Font for captions
        canvas.setFillColor('Grey')
        header_text = short_name+' DV Report, Compiled on '+date_str+' at '+time_str+' with PINES Analysis Toolkit V. '+str(pat_version)
        canvas.drawCentredString(612*0.5, 792-15, header_text)

    def page_footer():
        canvas.setFont('Times-Italic', 8) #Font for captions
        canvas.setFillColor('Grey')
        footer_text = str(page_num)
        canvas.drawCentredString(612*0.95, 792*0.02, footer_text)

    def title_page():
        page_header()
        global fig_num, page_num
        canvas.setFont('Times-Roman', 20) #Font for captions
        canvas.setFillColor('Black')
        title_text = 'PINES DV Report'
        canvas.drawCentredString(fig_x*0.5, fig_y*0.9, title_text)
        object_text = target+' ('+short_name+')'
        canvas.setFont('Times-Roman', 16) #Font for captions
        canvas.drawCentredString(fig_x*0.5, fig_y*0.87, object_text)
         #Add the source image. 
        source_image_path = source_path/('target_and_refs.png')
        img = ImageReader(source_image_path)
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

    def night_page(): 
        global fig_num, page_num
        #Get all nightly-normalized target lightcurve plots in the analysis directory.
        nightly_glob = natsorted(np.array([x for x in analysis_path.glob('*nightly*target_flux.png')]))

        if os.path.exists(analysis_path/('optimal_aperture.txt')):
            with open(analysis_path/('optimal_aperture.txt'), 'r') as f:
                lines = f.readlines()
                best_ap = lines[0].split('  ')[1].split('\n')[0]
                ap_type = best_ap.split('_')[1]
        best_aper_phot_path = analysis_path/('aper_phot_analysis/'+best_ap)

        nightly_glob = np.array(natsorted([x for x in best_aper_phot_path.glob('*nightly*.png')]))
        #Sort plots in order we want them in the report: arget lc, raw flux, normalized flux
        nightly_glob = np.array([nightly_glob[2], nightly_glob[1], nightly_glob[0]])

        if ap_type == 'fixed':
            rad = best_ap.split('_')[0]
            cap_ender = 'radius = {} pixels'.format(rad)
        elif ap_type == 'variable':
            fact = best_ap.split('_')[0]
            cap_ender = 'multiplicative seeing factor = {}'.format(fact)
            page_footer

        caption_texts = ['Best nightly-normalized corrected target flux, {} aperture photometry, {}.'.format(ap_type, cap_ender),
                        'Raw flux, {} aperture photometry, {}.'.format(ap_type, cap_ender),
                        'Nightly-normalized flux, {} aperture photometry, {}.'.format(ap_type, cap_ender)]
        
        for i in range(len(nightly_glob)):
            if i == 0:
                page_header()
            canvas.setFont('Times-Roman', 12) #Font for captions
            canvas.setFillColor('Black')
            image_path = str(nightly_glob[i])
            aperture_radius = image_path.split('=')[1].split('_')[0]
            
            img = ImageReader(image_path)
            y = (3-(i+1))*fig_y/3 + 50 
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

    def global_page(): 
        global fig_num, page_num

        if os.path.exists(analysis_path/('optimal_aperture.txt')):
            with open(analysis_path/('optimal_aperture.txt'), 'r') as f:
                lines = f.readlines()
                best_ap = lines[1].split(' ')[1]
                ap_type = best_ap.split('_')[1]
        best_aper_phot_path = analysis_path/('aper_phot_analysis/'+best_ap)
        

        global_glob = np.array(natsorted([x for x in best_aper_phot_path.glob('*global*.png')]))
        #Sort plots in order we want them in the report: arget lc, raw flux, normalized flux
        global_glob = np.array([global_glob[2], global_glob[1], global_glob[0]])

        if ap_type == 'fixed':
            rad = best_ap.split('_')[0]
            cap_ender = 'radius = {} pixels'.format(rad)
        elif ap_type == 'variable':
            fact = best_ap.split('_')[0]
            cap_ender = 'multiplicative seeing factor = {}'.format(fact)
            page_footer

        caption_texts = ['Best globally-normalized corrected target flux, {} aperture photometry, {}.'.format(ap_type, cap_ender),
                        'Raw flux, {} aperture photometry, {}.'.format(ap_type, cap_ender),
                        'Globally-normalized flux, {} aperture photometry, {}.'.format(ap_type, cap_ender)]
        
        for i in range(len(global_glob)):
            if i == 0:
                page_header()
            canvas.setFont('Times-Roman', 12) #Font for captions
            canvas.setFillColor('Black')
            image_path = str(global_glob[i])
            aperture_radius = image_path.split('=')[1].split('_')[0]
            
            img = ImageReader(image_path)
            y = (3-(i+1))*fig_y/3 + 50 
            w = fig_x #All the way across the page
            h = w * lc_aspect_ratio
            canvas.drawImage(img, x, y, x+w, h)
            caption_text = 'Figure {}: '.format(fig_num)+caption_texts[i]
            canvas.drawCentredString(x+w*0.5,y-h*0.1, caption_text)
            fig_num += 1

            if i == len(global_glob) - 1:
                page_footer()

        canvas.showPage()
        page_num += 1
    
    def centroid_page():
        global fig_num, page_num
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

    def diagnostics_page():
        global fig_num, page_num
        page_header()

        canvas.setFont('Times-Roman', 12) #Font for captions
        canvas.setFillColor('Black')
        seeing_path = diagnostic_plot_path/(short_name+'_seeing.png')
        img = ImageReader(seeing_path)
        y = 2*fig_y/3 +50
        w = fig_x #All the way across the page
        h = w * lc_aspect_ratio
        canvas.drawImage(img, x, y, x+w, h)
        caption_text = 'Figure {}: '.format(fig_num)+'Seeing measurements (in arcsec).'
        canvas.drawCentredString(x+w*0.5,y-h*0.1, caption_text)
        fig_num += 1


        background_path = diagnostic_plot_path/(short_name+'_backgrounds.png')
        img = ImageReader(background_path)
        y = 1*fig_y/3 + 50
        w = fig_x #All the way across the page
        h = w * lc_aspect_ratio
        canvas.drawImage(img, x, y, x+w, h)
        caption_text = 'Figure {}: '.format(fig_num)+'Background measurements (in ADU). Non-linear effects begin near 4000 ADU.'
        canvas.drawCentredString(x+w*0.5,y-h*0.1, caption_text)
        fig_num += 1

        airmass_path = diagnostic_plot_path/(short_name+'_airmasses.png')
        img = ImageReader(airmass_path)
        y = 0*fig_y/3 + 50
        w = fig_x #All the way across the page
        h = w * lc_aspect_ratio
        canvas.drawImage(img, x, y, x+w, h)
        caption_text = 'Figure {}: '.format(fig_num)+'Airmass measurements.'
        canvas.drawCentredString(x+w*0.5,y-h*0.1, caption_text)
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
        best_aper_phot_path = analysis_path/('aper_phot_analysis/'+best_ap+'/corr_ref_plots/')
        corr_glob = np.array(natsorted(list([x for x in best_aper_phot_path.glob('*.png')])))

        count = 0 
        for i in range(len(corr_glob)):
            if count % 3 == 0:
                page_header()
            canvas.setFont('Times-Roman', 12) #Font for captions
            canvas.setFillColor('Black')
            image_path = str(corr_glob[i])
            
            img = ImageReader(image_path)
            y = (3-((count%3)+1))*fig_y/3 + 50 
            w = fig_x #All the way across the page
            h = w * lc_aspect_ratio
            canvas.drawImage(img, x, y, x+w, h)
            if i == 0:
                caption_text = 'Figure {}: '.format(fig_num)+image_path.split('/')[-1].split('_')[0]+' corrected flux.'
            else:
                caption_text = 'Figure {}: '.format(fig_num)+'Reference '+str(i)+' corrected flux.'

            canvas.drawCentredString(x+w*0.5,y-h*0.1, caption_text)
            fig_num += 1
            
            count += 1
            if count % 3 == 0:
                page_footer()
                canvas.showPage()
                page_num += 1        

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
    global fig_num, page_num, fig_x, fig_y, source_aspect_ratio, lc_aspect_ratio, centroid_aspect_ratio, x
    fig_num = 1 #Initialize 
    page_num = 1 
    fig_x = 612
    fig_y = 792 #Dimensions (in pixels?) for 8.5x11 page
    source_aspect_ratio = 9/10
    lc_aspect_ratio = 5/17
    centroid_aspect_ratio = 27/51
    x = 0 #Start all figures at the left margin

    #------------------------------------------------------------------------------------------------------------------------
    #PAGE 1: Title and Source Detection Image.
    #Make the title page. 
    title_page()


    #------------------------------------------------------------------------------------------------------------------------
    #PAGE 2: Best nightly-normalized lightcurve, raw flux, and normalized flux. 
    night_page()
    
    #------------------------------------------------------------------------------------------------------------------------
    #PAGE 3: Best globally-normalized lightcurve, raw flux, and normalized flux. 
    global_page()

    #------------------------------------------------------------------------------------------------------------------------
    
    #PAGE 4: Centroid diagnostic plots. 
    centroid_page() 

    #------------------------------------------------------------------------------------------------------------------------
    #PAGE 5: Diagnostic plots.
    diagnostics_page()
    
    #------------------------------------------------------------------------------------------------------------------------
    #PAGES 6-N: Plots of corrected target flux. 
    corr_refs_pages()

    canvas.save()
    
    return

if __name__ == '__main__':
    target = '2MASS J00144919-0838207' 
    dv_report(target)