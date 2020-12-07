#!/usr/bin/env python

import setuptools
from numpy.distutils.core import Extension, setup

setup(
    name="VAISD",
    version="1.0",
    author="Daniele Padula",
    author_email="dpadula85@yahoo.it",
    description="A python package to compute intramolecular El-Ph coupling",
    url="https://github.com/dpadula85/VAISD",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU GPL License",
        "Operating System :: OS Independent",
    ],
	scripts=['VAISD/bin/compute_specden'],
    zip_safe=False
)
