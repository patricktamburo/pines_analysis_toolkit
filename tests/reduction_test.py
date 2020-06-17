import pines_analysis_toolkit as pat

#List of 202003 target names
targets = ['2MASS J11250502+0556424', '2MASS J11394192-0310039', '2MASS J11491231-0153003', '2MASS J11555389+0559577', '2MASS J11582077+0435014']

for target in targets:
    pat.data.get_raw_science_files(target)

    pat.data.reduce(target, upload=True, delete_raw=True, delete_reduced=False)