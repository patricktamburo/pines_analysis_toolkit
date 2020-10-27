import numpy as np
'''Authors:
		Patrick Tamburo, Boston University, June 2020
	Purpose:
        Identifies the id of the target given an initial guess of its location. 
	Inputs:
        sources (pandas dataframe): Table of detected sources from detect_sources.
        guess_position (tuple of floats): Guess of targets position (x pixel, y pixel) in the image for which source detection was performed.
    Outputs:
        target_id (int): the id of the target in the source table that most closely matches guess_position
	TODO:
'''
def target_finder(sources, guess_position):
    distances = np.sqrt((guess_position[0] - sources['xcenter'])**2 + (guess_position[1] - sources['ycenter'])**2)
    target_id = np.where(distances == np.min(distances))[0][0]
    return target_id