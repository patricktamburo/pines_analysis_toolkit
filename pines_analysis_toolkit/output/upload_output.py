from numpy.lib.utils import source
from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
from pines_analysis_toolkit.utils.pines_login import pines_login
import os 
from natsort import natsorted
import numpy as np 
from pathlib import Path

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

    #Get the path of the object's DV report on the local machine. 
    dv_report_path = pines_path/('Objects/'+short_name+'/output/'+short_name.replace(' ','')+'_dv_report.pdf')

    #Move the sftp connection to the object's lightcurves/ directory on the PINES server. 
    sftp.chdir('/data/lightcurves/')
    #If the folder doesn't exist, make it. 
    if not short_name in sftp.listdir():
        sftp.mkdir(short_name)
    sftp.chdir(short_name)

    #Upload the DV report
    if os.path.exists(dv_report_path):
        print('Uploading {} DV report.'.format(short_name))
        upload_path = sftp.pwd+'/'+short_name.replace(' ','')+'_dv_report.pdf'
        sftp.put(dv_report_path, upload_path)
    else:
        print('WARNING: No DV report found for {}.'.format(short_name))

    if not 'sources' in sftp.listdir():
        sftp.mkdir('sources')

    sftp.chdir('sources')

    #Upload source detection CSV.
    source_detection_path = pines_path/('Objects/'+short_name+'/sources/target_and_references_source_detection.csv')
    if os.path.exists(source_detection_path):
        print('Uploading {} source detection CSV.'.format(short_name))
        upload_path = sftp.pwd+'/target_and_references_source_detection.csv'
        sftp.put(source_detection_path, upload_path)
    else:
        print('WARNING: No source detection CSV found for {}.'.format(short_name))
    
    #Upload centroids CSV. 
    centroid_path = pines_path/('Objects/'+short_name+'/sources/target_and_references_centroids.csv')
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
    aper_phot_files = np.array(natsorted([x for x in aper_phot_path.glob('*.csv')]))
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

if __name__ == '__main__':
    sftp = pines_login()
    targets = ['2MASS J09161504+2139512',
'2MASS J09183815+2134058',
'2MASS J09202645+2651583',
'2MASS J09230861+2340152',
'2MASS J09153413+0422045',
'2MASS J09165865+1057192',
'2MASS J09214255+0842030',
'2MASS J09385888+0443438',
'2MASS J09413492+1009428',
'2MASS J06420559+4101599',
'2MASS J11250502+0556424',
'2MASS J11394192-0310039',
'2MASS J11491231-0153003',
'2MASS J11555389+0559577',
'2MASS J11582077+0435014',
'2MASS J12403680+1525172',
'2MASS J13004255+1912354',
'2MASS J13043318+0907070',
'2MASS J13054106+2046394',
'2MASS J13083106+0818522',
'2MASS J13154557+2454072',
'2MASS J13340623+1940351',
'2MASS J12321827-0951502',
'2MASS J12545393-0122474',
'2MASS J12565688+0146163',
'2MASS J12573726-0113360',
'2MASS J13334540-0215599',
'2MASS J14225717+0827539',
'2MASS J14243909+0917104',
'2MASS J14313097+1436539',
'2MASS J14391190+0823157',
'2MASS J14392836+1929149',
'2MASS J14482563+1031590',
'2MASS J14520183+1114590',
'2MASS J15065441+1321060',
'2MASS J16202614-0416315',
'2MASS J16235981-0508066',
'2MASS J16291840+0335371',
'2MASS J16304139+0938446',
'2MASS J16360078-0034525',
'2MASS J17312974+2721233',
'2MASS J17333281+3144571',
'2MASS J17392515+2454421',
'2MASS J17434148+2127069',
'2MASS J17561080+2815238',
'2MASS J17260007+1538190',
'2MASS J17281134+0839590',
'2MASS J17433635+1549025',
'2MASS J18212815+1414010',
'2MASS J20304235+0749358',
'2MASS J20343769+0827009',
'2MASS J20360316+1051295',
'2MASS J20575409-0252302',
'2MASS J21073169-0307337',
'2MASS J21304464-0845205',
'2MASS J21373742+0808463',
'2MASS J21392676+0220226',
'2MASS J22081363+2921215',
'2MASS J22425317+2542573',
'2MASS J22443167+2043433',
'2MASS J22490917+3205489',
'2MASS J22244381-0158521',
'2MASS J22552907-0034336',
'2MASS J23062928-0502285',
'2MASS J23255604-0259508',
'2MASS J00144919-0838207',
'2MASS J00191165+0030176',
'2MASS J00242463-0158201',
'2MASS J00261147-0943406',
'2MASS J00320509+0219017',
'2MASS J00345684-0706013',
'2MASS J00464841+0715177',
'2MASS J00492826+0440575',
'2MASS J01075242+0041563',
'2MASS J01221697+0331235',
'2MASS J01353586+1205216',
'2MASS J01365662+0933473',
'2MASS J01410321+1804502',
'2MASS J01535422+1404532',
'2MASS J02050344+1251422',
'2MASS J02073557+1355564',
'2MASS J04132039-0114248',
'2MASS J04223057+0723448',
'2MASS J04234858-0414035',
'2MASS J04433761+0002051',
'2MASS J04574602-0207179',
'2MASS J05002100+0330501',
'2MASS J05012406-0010452',
'2MASS J08312221+1538511',
'2MASS J08350622+1953050',
'2MASS J08354537+2224310',
'2MASS J08433323+1024470',
'2MASS J08560211+1240150',
'2MASS J10151763+0752121',
'2MASS J10292165+1626526',
'2MASS J10330903+1216265',
'2MASS J10433508+1213149',
'2MASS J10515129+1311164',
'2MASS J15394189-0520428',
'2MASS J15551573-0956055',
'2MASS J15552614+0017204',
'2MASS J16051015+0421521',
'2MASS J16192830+0050118']
    for target in targets:
        upload_output(sftp, target)