import pines_analysis_toolkit as pat

#List of 202004 target names
targets = ['2MASS J13054106+2046394','2MASS J13154557+2454072',
           '2MASS J13340623+1940351','2MASS J13083106+0818522','2MASS J12321827-0951502','2MASS J12545393-0122474','2MASS J12565688+0146163',
           '2MASS J12573726-0113360','2MASS J13334540-0215599']

for target in targets:
    pat.data.get_raw_science_files(target)

    pat.data.reduce(target, upload=True, delete_raw=True, delete_reduced=False)