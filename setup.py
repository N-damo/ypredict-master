from setuptools import setup
from setuptools import find_packages
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'),'r') as f:
    long_description = f.read()

setup(name='Ypredict',version='0.1.3',description='a very fast Y-haplogroup caller',long_description=long_description, url='https://github.com/N-damo/ypredict-master',
author="Li'an Lin",author_email='21620151153308@stu.xmu.edu.cn',license='GPL3', keywords='Y haplogroup', # have to be included in MANIFEST.in as well.
package_data={'test':['./test/hfspecial.xlsx','./test/vcfhead','./test/map.json','./test/snp14.3.csv','./test/ref_vcf.gz','./test/filter.csv']}
,include_package_data=True, packages=['bin'],install_requires=['biopython >=1.72'],scripts=['bin/ypredict.py','bin/snpfilter.py'])
