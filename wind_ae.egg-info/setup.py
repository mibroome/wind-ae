#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup

setup(
    name="wind_ae",
    version="1.0",
    author='Madelyn Broome',
    author_email='mabroome@ucsc.edu',
    description="1D relaxation Parker wind model based on Murray-Clay et al. (2009)",
    install_requires=[line.strip() for line in
                      open('requirements.txt', 'r').readlines()],
    packages=find_packages(),
    package_data={
        'wind_ae': ['wind_ae/wrapper/wrapper_utils/*.dat',
                    'wind_ae/saves/*.csv',
                    'wind_ae/inputs/*', 
                    'wind_ae/McAstro/atoms/*.csv',
                    'wind_ae/McAstro/atoms/*.txt',
                    'wind_ae/McAstro/atoms/Lodders.dat'],  # include all .txt files in the data folder
    },
    include_package_data=True,
)