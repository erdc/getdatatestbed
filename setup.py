#!/usr/bin/env python
from setuptools import setup

setup(name='GetDataTestbed',
      version='0.1.dev',
      description=('Utilities for retrieving data relevant to the USACE' +
                   'Field Research Facility Coastal Model Test Bed' +
                   '(CMTB)'),
      author='Spicer Bak',
      modules=['getDataFRF', 'getOutsideData', 'download_grid_data'],
     )