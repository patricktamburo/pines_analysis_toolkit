import pandas as pd 
import pdb 
import distutils

def profile_reader(profile_path):
    df = pd.read_csv(profile_path, sep='=', header=None).transpose()
    new_header = df.iloc[0]
    df = df[1:]
    df.columns = new_header
    df.columns = df.columns.str.strip()

    output_dict = {'source_detect_image':df['source_detect_image'].iloc[0].strip(),
                   'exclude_lower_left':bool(distutils.util.strtobool(df['exclude_lower_left'].iloc[0].strip())),
                   'dimness_tolerance':float(df['dimness_tolerance'].iloc[0]),
                   'brightness_tolerance':float(df['brightness_tolerance'].iloc[0]),
                   'distance_from_target':float(df['distance_from_target'].iloc[0]),
                   'non_linear_limit':float(df['non_linear_limit'].iloc[0]),
                   'edge_tolerance':float(df['edge_tolerance'].iloc[0]),
                   'guess_position_x':float(df['guess_position_x'].iloc[0]),
                   'guess_position_y':float(df['guess_position_y'].iloc[0])}
                   
    return output_dict