# -*- coding: utf-8 -*-

from distutils.core import setup
import os
import shutil
shutil.copy('README.md', 'tisutil/README.md')

dir_setup = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(dir_setup, 'tisutil', 'release.py')) as f:
    # Defines __version__
    exec(f.read())

setup(name='tisutil',
      version=__version__,
      description='Extension of matfuncutil for scattering matrices and other quantities.',
      author="Peter Bingham",
      author_email="petersbingham@hotmail.co.uk",
      packages=['tisutil'],
      package_data={'tisutil': ['tests/*', 'README.md']}
     )
