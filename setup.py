#!/usr/bin/env python

from distutils.core import setup

setup(name='GetDataTestbed',
      version='0.1dev',
      description=('Utilities for retrieving data relevant to the USACE' +
                   'Field Research Facility Computational Model Test Bed' +
                   '(CMTB)'),
      author='Spicer Bak',
      modules=['getDataFRF', 'getOutsideData', 
               'download_grid_data'],
     )
