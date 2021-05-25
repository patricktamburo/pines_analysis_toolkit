import pines_analysis_toolkit as pat
import pdb 
from pathlib import Path

target = '2MASS 23060502_3'
pines_path = pat.utils.pines_dir_check()
short_name = pat.utils.short_name_creator(target)

if target == '2MASS 23060502_2':
    manual_flat_path = Path('/Users/tamburo/Documents/PINES_analysis_toolkit/Calibrations/Flats/Domeflats/J/Master Flats/master_flat_J_20201003.fits')
elif target == '2MASS 23060502_3':
    manual_flat_path = Path('/Users/tamburo/Documents/PINES_analysis_toolkit/Calibrations/Flats/Domeflats/J/Master Flats/master_flat_J_20200810.fits')
else:
    raise ValueError('target must be "2MASS 23060502_2" or "2MASS 23060502_3"')

pat.data.reduce(target, manual_flat_path=manual_flat_path)
profile_data = pat.utils.profile_reader(pines_path/('Objects/'+short_name+'/2mass2306-0502.profile'))
sources = pat.photometry.ref_star_chooser(target, profile_data['source_detect_image'], restore=True)
centroided_sources = pat.photometry.centroider(target, sources, restore=True)
pat.photometry.aper_phot.fixed_aper_phot(target, centroided_sources, [4.0, 4.5])

pat.analysis.weighted_lightcurve(target, plots=True)
pdb.set_trace()