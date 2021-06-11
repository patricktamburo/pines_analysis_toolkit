import pines_analysis_toolkit as pat 

targets = ['2MASS J15394189-0520428'
'2MASS J15551573-0956055',
'2MASS J15552614+0017204',
'2MASS J16051015+0421521',
'2MASS J16192830+0050118']

sftp = pat.utils.pines_login()
for target in targets:
    pat.data.get_raw_science_files(sftp, target)
    pat.data.reduce(target, delete_raw=False, delete_reduced=False, linearity_correction=True, linearity_correction_degree=4)

sftp.close()

 