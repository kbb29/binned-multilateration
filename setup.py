#!/usr/bin/env python

from distutils.core import setup

setup(name='binned_multilateration',
      version='1.0',
      description='Multilateration for binned distances',
      author='kbb29',
      author_email='not@nemail.addr',
      url='https://github.com/kbb29/bucket-multilateration',
      install_requires=['scipy', 'numpy', 'chevron'],
      packages=['binned_multilateration'],
      package_data={'binned_multilateration': ['map.mustache'] }
     )
