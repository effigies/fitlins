# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Base module variables
"""

__version__ = '0.0.1-dev'
__author__ = 'The CRN developers'
__copyright__ = 'Copyright 2017, Center for Reproducible Neuroscience, Stanford University'
__credits__ = ['Christopher J. Markiewicz', 'Chris Gorgolewski', 'Russell A. Poldrack']
__license__ = '3-clause BSD'
__maintainer__ = 'Christopher J. Markiewicz'
__email__ = 'crn.poldracklab@gmail.com'
__status__ = 'Prototype'
__url__ = 'https://github.com/poldracklab/fitlins'
__packagename__ = 'fitlins'
__description__ = ''
__longdesc__ = ''

DOWNLOAD_URL = (
    'https://github.com/effigies/{name}/archive/{ver}.tar.gz'.format(
        name=__packagename__, ver=__version__))


SETUP_REQUIRES = [
    'setuptools>=18.0',
]

REQUIRES = [
    'nistats',
    'pybids>=0.3',
    'niworkflows>=0.1.8',
]

LINKS_REQUIRES = [
]

TESTS_REQUIRES = [
]

EXTRA_REQUIRES = {
    'doc': ['sphinx>=1.5.3', 'pydotplus', 'pydot>=1.2.3', 'sphinx_rtd_theme', 'sphinx-argparse'],
    'tests': TESTS_REQUIRES,
    'duecredit': ['duecredit'],
}

# Enable a handle to install all extra dependencies at once
EXTRA_REQUIRES['all'] = [val for _, val in list(EXTRA_REQUIRES.items())]
CLASSIFIERS = [
    'Development Status :: 2 - Pre-Alpha',
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering',
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
]
