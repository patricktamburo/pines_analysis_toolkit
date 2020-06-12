from pines_analysis_toolkit.utils.short_name_creator import short_name_creator
from pines_analysis_toolkit.utils.pines_dir_check import pines_dir_check
import pdb

def detect_sources(target):
    #Get the target's 'short name'
    short_name = short_name_creator(target)

    #Get your local PINES directory
    pines_path = pines_dir_check()

    pdb.set_trace()
