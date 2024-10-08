from setuptools import setup
from setuptools import find_packages

setup( name = 'pyStrint',
version = '0.0.13',
description='Deciphering more accurate cell-cell interactions by modeling cells and their interactions',
url='https://github.com/deepomicslab/StrInt',
author='Jingwan WANG',
author_email='wanwang6-c@my.cityu.edu.hk',
license='MIT',
packages=find_packages(),
install_requires = [
    'numpy==1.22.3',
    'umap-learn==0.5.2',
    'loess==2.1.2',
    'pandas==1.5.2',
    'scipy==1.9.3',
    'scanpy==1.9.8',
    'smurf-imputation'
],
package_data={'LR': ['*.txt'],'pipelines':['*']},
include_package_data=True
)