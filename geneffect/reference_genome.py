from __future__ import absolute_import, division, print_function

import os

from .util import as_biopython_seq

class GenomeReader(object):
    
    def __init__(self, config_setup):
        self.ref_genome_dir = config_setup.get_path('REFERENCE_GENOME_DIR')
        self.chromosome_readers = dict(map(self._create_chromosome_reader, os.listdir(self.ref_genome_dir)))

    def read_seq(self, chromosome, start, end):
        return self.chromosome_readers[find_chrom(chromosome, self.chromosome_readers.keys())].read_seq(start, end)
        
    def close(self):
        for chromosome_reader in self.chromosome_readers.values():
            chromosome_reader.file_handler.close()
            
    def __contains__(self, chromosome):
        return find_chrom(chromosome, self.chromosome_readers.keys()) is not None
        
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

def find_chrom(query_chr_name, available_chr_names):

    assert isinstance(query_chr_name, str), 'Unexpected chromosome type: %s' % type(query_chr_name)

    if query_chr_name.lower().startswith('chr'):
        query_chr_name = query_chr_name[3:]
        
    query_chr_name = query_chr_name.upper()
    
    for possible_chr_name in _find_synonymous_chr_names(query_chr_name):
        
        if possible_chr_name in available_chr_names:
            return possible_chr_name
            
        prefixed_possible_chr_name = 'chr%s' % possible_chr_name
            
        if prefixed_possible_chr_name in available_chr_names:
            return prefixed_possible_chr_name
            
    return None
        
def _find_synonymous_chr_names(chr_name):
    
    for synonymous_chr_name_group in _SYNONYMOUS_CHR_NAME_GROUPS:
        if chr_name in synonymous_chr_name_group:
            return synonymous_chr_name_group
            
    # Single-digit numbers can either appear with or without a trailing 0.
    if chr_name.isdigit() and len(str(int(chr_name))) == 1:
        chr_number = str(int(chr_name))
        return {chr_number, '0' + chr_number}
            
    return {chr_name}

_SYNONYMOUS_CHR_NAME_GROUPS = [
    {'X', '23'},
    {'Y', '24'},
    {'XY', '25'},
    {'M', 'MT', '26'},
]
