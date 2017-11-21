from __future__ import absolute_import, division, print_function

from .util import log
from .config import ConfigSetup
from .reference_genome import GenomeReader
from .uniprot_and_pfam_data_loader import load_pfam_data, load_uniprot_data
from .gene_data_loader import load_genes
from .variant_processing import VariantInterpreter

class Setup(object):
    def __init__(self, user_specified_ref_genome, should_load_pfam_data = True):
        
        self._config_setup = ConfigSetup(user_specified_ref_genome)
        self.genome_reader = GenomeReader(self._config_setup)
        
        if should_load_pfam_data:
            self.pfam_data = load_pfam_data(self._config_setup)
        else:
            log('Skipping on setting pfam data.')
            self.pfam_data = None
            
        self.uniprot_records = load_uniprot_data(self._config_setup, self.pfam_data)
        self.genes = load_genes(self._config_setup, self.genome_reader, self.uniprot_records)
        self.variant_interpreter = VariantInterpreter(self.genome_reader, self.genes)
        
        log('Finished setting up reference genome %s.' % (self._config_setup.ref_genome))
