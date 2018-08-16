# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 20:54:36 2018

@author: Martin
"""

from setuptools import setup, find_packages
import os

MAJOR = 0
MINOR = 0
MICRO = 1
ISRELEASED = False
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

DESCRIPTION = ("Tools for performing kinetic modelling of biomolecular networks")

def write_version_py(filename=r'py_digilab/version.py'):

    cnt = """
# THIS FILE IS GENERATED FROM SCIPY SETUP.PY
short_version = '{version}'
version = '{version}'
release = {isrelease}
"""
    path = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(path,filename), 'w+') as a:
        a.write(cnt.format(**{'version': VERSION,
                       'isrelease': str(ISRELEASED)}))

write_version_py()

setuptools_kwargs = {
    'name': 'py-digilab',
    'version': VERSION,
    'author': 'Martin Wong',
    'description': DESCRIPTION,
    'packages':find_packages(exclude=['tests.*','tests*']),
    'install_requires': [
        'scipy',
        'numpy',
        'pandas'
    ],
    'setup_requires': ['numpy'],
    'zip_safe': False,
}

setup(**setuptools_kwargs)