from __future__ import absolute_import, division, print_function

import os

from .util import as_biopython_seq

class GenomeReader(object):
    
    def __init__(self, config_setup):
        self.ref_genome_dir = config_setup.get_path('REFERENCE_GENOME_DIR')
        self.chromosome_readers = dict(map(self._create_chromosome_reader, os.listdir(self.ref_genome_dir)))

    def read_seq(self, chromosome, start, end):
        return self.chromosome_readers[_fix_chromosome_name(chromosome)].read_seq(start, end)
        
    def close(self):
        for chromosome_reader in self.chromosome_readers.values():
            chromosome_reader.file_handler.close()
            
    def __contains__(self, chromosome):
        return _fix_chromosome_name(chromosome) in self.chromosome_readers
        
    def _create_chromosome_reader(self, file_name):
        chr_name = file_name.split('.')[0].replace('chr', '')
        f = open(os.path.join(self.ref_genome_dir, file_name), 'r')
        return (chr_name, ChromosomeReader(f))

class ChromosomeReader(object):

    def __init__(self, file_handler):
        self.file_handler = file_handler
        self.header_len = len(file_handler.readline())
        self.line_len = len(file_handler.readline()) - 1

    def read_seq(self, start, end):

        absolute_start = self.convert_to_absolute_coordinate(start)
        absolute_length = self.convert_to_absolute_coordinate(end) - absolute_start + 1

        self.file_handler.seek(absolute_start)
        seq = self.file_handler.read(absolute_length).replace('\n', '').upper()
        return as_biopython_seq(seq)
        
    def convert_to_absolute_coordinate(self, position):
        position_zero_index = position - 1
        return self.header_len + position_zero_index + (position_zero_index // self.line_len)
    
def _fix_chromosome_name(chromosome):
    if chromosome == 'MT':
        return 'M'
    else:
        return chromosome
