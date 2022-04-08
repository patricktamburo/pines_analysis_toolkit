import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='pines_analysis_toolkit',
    version='1.0.0',
    author='Patrick Tamburo',
    author_email='tamburop@bu.edu',
    description='Testing installation of Package',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/patricktamburo/pines_analysis_toolkit',
    license='MIT',
    packages=['pines_analysis_toolkit'],
    install_requires=['julian'],
)