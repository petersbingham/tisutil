# -*- coding: utf-8 -*-

from distutils.core import setup
import shutil
shutil.copy('README.md', 'tisutil/README.md')

setup(name='tisutil',
      version='0.6',
      description='Extension of matfuncutil for scattering matrices and other quantities.',
      author="Peter Bingham",
      author_email="petersbingham@hotmail.co.uk",
      packages=['tisutil'],
      package_data={'tisutil': ['tests/*', 'README.md']}
     )
