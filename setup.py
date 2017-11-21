from __future__ import absolute_import, division, print_function

import os
from shutil import copyfile

from setuptools import setup
from setuptools.command.install import install

DEFAULT_CONFIG_FILE_PATH = 'geneffect/default_config.py'
OPERATIONAL_CONFIG_FILE_PATH = os.path.expanduser('~/.geneffect_config.py')

class CreateConfigFileCommand(install):
    def run(self):
        copyfile(DEFAULT_CONFIG_FILE_PATH, OPERATIONAL_CONFIG_FILE_PATH)
        print('Created your configuration file with default setting at: %s. Please review those settings.' % OPERATIONAL_CONFIG_FILE_PATH)

setup(
    name = 'geneffect',
    version = '1.0',
    description = 'A Python library for retrieving functional annotations of genes and analyzing the ' + \
            'effects of genetic variants, currently focusing on proteomic data of protein-coding genes.',
    url = 'https://github.com/nadavbra/geneffect',
    author = 'Nadav Brandes',
    author_email  ='nadav.brandes@mail.huji.ac.il',
    license = 'MIT',
    packages = ['geneffect'],
    install_requires = [
        'numpy',
        'pandas',
        'biopython',
        'interval_tree',
    ],
    cmdclass = {
        'install': CreateConfigFileCommand,
    },
)
