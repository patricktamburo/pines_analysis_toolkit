import os
from pathlib import Path
import pdb
import setuptools
import time

with open('README.md', 'r') as fh:
    long_description = fh.read()

home_dir = Path(os.path.expanduser('~/Documents/'))
default_pines_dir = os.path.join(home_dir,'PINES_analysis_toolkit/')
print('')
print('')
print('')
print('Setting up analysis directory at:',default_pines_dir)
print('All PINES data products will be stored in this directory.')
time.sleep(5)

if not os.path.exists(default_pines_dir):
    print('Creating PINES directory at',str(default_pines_dir))
    os.mkdir(default_pines_dir)
    os.mkdir(os.path.join(default_pines_dir,'Calibrations/'))
    os.mkdir(os.path.join(default_pines_dir,'Calibrations/Flats/'))
    os.mkdir(os.path.join(default_pines_dir,'Calibrations/Flats/Domeflats/'))
    os.mkdir(os.path.join(default_pines_dir,'Calibrations/Flats/Domeflats/J/'))
    os.mkdir(os.path.join(default_pines_dir,'Calibrations/Flats/Domeflats/J/Raw/'))
    os.mkdir(os.path.join(default_pines_dir,'Calibrations/Flats/Domeflats/J/Master Flats/'))
    os.mkdir(os.path.join(default_pines_dir,'Calibrations/Flats/Domeflats/H/'))
    os.mkdir(os.path.join(default_pines_dir,'Calibrations/Flats/Domeflats/H/Raw/'))
    os.mkdir(os.path.join(default_pines_dir,'Calibrations/Flats/Domeflats/H/Master Flats/'))
    os.mkdir(os.path.join(default_pines_dir,'Calibrations/Darks/'))
    os.mkdir(os.path.join(default_pines_dir,'Calibrations/Darks/Raw/'))
    os.mkdir(os.path.join(default_pines_dir,'Calibrations/Darks/Master Darks/'))
    os.mkdir(os.path.join(default_pines_dir,'Calibrations/Bad Pixel Masks/'))
    os.mkdir(os.path.join(default_pines_dir,'Logs/'))
    os.mkdir(os.path.join(default_pines_dir,'Objects/'))
else:
    print(default_pines_dir,'already exists, skipping.')

setuptools.setup(
    name='pines_analysis_toolkit', 
    version='0.0.1',
    author='Patrick Tamburo',
    author_email='tamburop@bu.edu',
    description='A package for accessing, reducing, and analyzing PINES data.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='',
    packages=setuptools.find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
)



