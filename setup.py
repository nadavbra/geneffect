from setuptools import setup

setup(name = 'geneffect',
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
      ])
