from setuptools import setup
from setuptools import find_packages

setup(name='Ypredict',version='0.0.1',packages=['bin'],install_requires=['biopython >=1.72'],scripts=['bin/ypredict.py','bin/snpfilter.py'])
