from setuptools import setup
import re


# borrowed from Agile Scientific bruges tools
verstr = 'unknown'
VERSIONFILE = "psm/_version.py"
with open(VERSIONFILE, "r") as f:
    verstrline = f.read().strip()
    pattern = re.compile(r"__version__ = ['\"](.*)['\"]")
    mo = pattern.search(verstrline)
if mo:
    verstr = mo.group(1)
    print("Version "+verstr)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))



setup(
    name='psm',
    version=verstr,
    author='Sylvia Dee and Amir Allam',
    author_email='...',
          packages=['psm',
          'psm.agemodels',
          'psm.coral',
          'psm.cellulose',
          'psm.speleo',
          'psm.icecore',
          'psm.aux_functions'],,
    url='https://github.com/sylvia-dee/PRYSM',
    license='LICENSE',
    description='Climate Proxy System Modeling Tools in Python',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    install_requires=[
        'scipy',
        'numpy'],
    include_package_data=True,
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Intended Audience :: Science/Research'],
    project_urls={
        'Bug Reports': 'https://github.com/sylvia-dee/PRYSM/issues',
        'Source': 'https://github.com/sylvia-dee/PRYSM',
    },
)
