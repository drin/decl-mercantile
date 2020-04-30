# -*- coding: utf-8 -*-
from setuptools import setup

packages     = ['skyhookdm_singlecell', 'skyhookdm_singlecell.Tables']
package_data = {'': ['*']}

install_requires = [
    'awscli          >= 1.17.15 , <  2.0.0',
    'cython          >= 0.29.15 , < 0.30.0',
    'flatbuffers     >= 1.11    , <   1.12',
    'h5py            >= 2.10    , <    3.0',
    'numpy           >= 1.17    , <    2.0',
    'pandas          >= 1.0.1   , <  2.0.0',
    'pyarrow         >= 0.17.0  , < 0.18.0',
    'python-cephlibs >= 0.94.5  , < 0.95.0',
    'scipy           >= 1.3     , <    2.0',
]

setup_kwargs = {
    'name'            : 'skyhookdm-singlecell',
    'version'         : '0.1.0',

    'description'     : '',
    'long_description': None,

    'author'          : 'Aldrin Montana',
    'author_email'    : 'akmontan@ucsc.edu',

    'maintainer'      : None,
    'maintainer_email': None,

    'url'             : None,

    'packages'        : packages,
    'package_data'    : package_data,

    'install_requires': install_requires,
    'python_requires' : '>=3.6,<4.0',
}


setup(**setup_kwargs)
