import os 
def object_directory_creator(pines_path, short_name):
    '''Sets up object directory layout in pat directory.'''
    os.mkdir(pines_path/('Objects/'+short_name))
    os.mkdir(pines_path/('Objects/'+short_name+'/analysis'))
    os.mkdir(pines_path/('Objects/'+short_name+'/analysis/diagnostic_plots'))
    os.mkdir(pines_path/('Objects/'+short_name+'/analysis/aper_phot_analysis'))
    os.mkdir(pines_path/('Objects/'+short_name+'/aper_phot'))
    os.mkdir(pines_path/('Objects/'+short_name+'/output'))
    os.mkdir(pines_path/('Objects/'+short_name+'/psf_phot'))
    os.mkdir(pines_path/('Objects/'+short_name+'/pwv'))
    os.mkdir(pines_path/('Objects/'+short_name+'/pwv'))
    os.mkdir(pines_path/('Objects/'+short_name+'/pwv/CONUS'))
    os.mkdir(pines_path/('Objects/'+short_name+'/sources'))
    os.mkdir(pines_path/('Objects/'+short_name+'/raw'))
    os.mkdir(pines_path/('Objects/'+short_name+'/reduced'))
    return