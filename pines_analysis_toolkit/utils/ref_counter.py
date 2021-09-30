import pines_analysis_toolkit as pat 
import pdb 
from natsort import natsorted
import os 
import pandas as pd 
import numpy as np 
from scipy.stats import sigmaclip

def ref_counter():
    pines_path = pat.utils.pines_dir_check()
    objects_path = pines_path/('Objects/')
    targets = natsorted(os.listdir(objects_path))
    targets = np.array(targets[0:-1])
    targets = np.delete(targets, np.where(targets == '2MASS 0438+4349')[0][0]) #Remove staring sdM target.
    targets = np.delete(targets, np.where(targets == '2MASS 1333-0215')[0][0]) #Remove target that is saturated.
    num_refs = np.zeros(len(targets))
    for i in range(len(targets)):
        source_path = objects_path/(targets[i]+'/sources/target_and_references_source_detection.csv')
        try:
            df = pd.read_csv(source_path)
            num_refs[i] = len(df) - 1 #Subtract the target. 
        except:
            num_refs[i] = np.nan
    
    num_refs = num_refs[~np.isnan(num_refs)]
    v,l,h = sigmaclip(num_refs, 3, 3)
    breakpoint()
    return 

if __name__ == '__main__':
    ref_counter()