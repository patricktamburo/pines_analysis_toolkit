import pines_analysis_toolkit as pat

#Declare a target.
target = '2MASS J16291840+0335371'

#Detect sources in the target's field. 
sources = pat.photometry.detect_sources(target)
