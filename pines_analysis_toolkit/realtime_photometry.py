import PINES as p
import pdb
from photometry_main import photometry_main

#Group 26
#target = '2MASSJ09130320+1841501' 
#target = '2MASSJ09161504+2139512' #Brightest one this run
#target = '2MASSJ09183815+2134058'
#target = '2MASSJ09202645+2651583'
#target = '2MASSJ09230861+2340152'

#Group 27
#target = '2MASSJ09153413+0422045'
target = '2MASSJ09165865+1057192'
#target = '2MASSJ09214255+0842030'
#target = '2MASSJ09385888+0443438'
#target = '2MASSJ09413492+1009428'

#target = '2MASSJ06420559+4101599' #Jackie and Johanna's variable Spitzer target

#Group 39 
#target = '2MASSJ11250502+0556424'
#target = '2MASSJ11394192-0310039'
#target = '2MASSJ11491231-0153003'
#target = '2MASSJ11555389+0559577'
#target = '2MASSJ11582077+0435014'


#Group 62
target = '2MASSJ16291840+0335371'

split_name = target.split('J')[1]
short_name = '2MASS '+split_name[0:4]+split_name[8]+split_name[9:13]
#local_path = 'P:\\Photometry\\PINES\\'+run_directory+'\\'+short_name+'\\'
p.data.move_science_files(target)
p.reduction.reduce(short_name,flat_type='dome')
photometry_main(short_name,[5,6,7],run_identify_sources=0 ,guess_position=(705,386))
p.analysis.lightcurve(short_name,'5')
pdb.set_trace()