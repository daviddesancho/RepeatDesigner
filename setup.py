#!/usr/bin/env python

# Setup script for RepeatDesign package

import os
from setuptools import setup,find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
		name='RepeatDesigner',
		version='0.1dev',
		description='A library for building protein repeats',
		url='http://github.com/daviddesancho/repeatdesigner',
		author='David De Sancho',
		author_email='daviddesancho.at.gmail.com',
		license='GNU Lesser General Public License',
		packages=find_packages(),
                install_requires=[
                    'numpy', 'matplotlib', 'biopython', 'seaborn', 'cython',
                    ],
		keywords= " protein repeat",
		long_description=read('README.md'),
		classifiers = ["""\
				Development Status :: 1 - Planning
				Operating System :: POSIX :: Linux
				Operating System :: MacOS
				Programming Language :: Python :: 2.7
				Topic :: Scientific/Engineering :: Bio-Informatics
				Topic :: Scientific/Engineering :: Chemistry
				"""]
		)
