#  ___________________________________________________________________________
#
# MPC4CLR: Model Predictive Control for Critical Load Restoration 
#          in Power Distribution Systems 
# Copyright 2022 National Renewable Energy Laboratory (NREL) 
# This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________


import glob
import sys
import os
from setuptools import find_packages
from distutils.core import setup

# Raise an error if trying to install with python2
if sys.version[0] == '2':
    print("Error: This package must be installed with python3")
    sys.exit(1)

setup(
    name='mpc4clr',
    version='0.0.1',
    description="MPC for Critical Load Restoration in Power Distribution Systems.",
    url='https://github.nrel.gov/AGM-CSSO/Optimization_for_Grid_Resiliency/tree/MPC4CLR', # TODO: Update this url when the repo moves to github.com
    author='Abinet Tesfaye Eseye, Xiangyu Zhang, Bernard Knueven, Matthew Reynolds, Weijia Liu, and Wesley Jones',
    maintainer_email='Wesley.Jones@nrel.gov', # TODO: update the email if needed
    license='Revised BSD',
    packages=find_packages(),  
    scripts=[],
    include_package_data=True,
    install_requires=['matplotlib', 'pandas', 'numpy', 'jupyter', 'notebook', 'pytest'],    
    python_requires='>=3.8'
    )