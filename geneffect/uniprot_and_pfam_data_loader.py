from __future__ import absolute_import, division, print_function

import re
import gzip

import pandas as pd

from Bio import SeqIO

from .util import log

MIN_COUNTS_FOR_COMMON_CLAN = 25

class PfamData(object):
    def __init__(self, pfam_data, pfam_domains, pfam_domains_by_uniprot_id, pfam_clan_counts, common_pfam_clans):
        self.pfam_data = pfam_data
        self.pfam_domains = pfam_domains
        self.pfam_domains_by_uniprot_id = pfam_domains_by_uniprot_id
        self.pfam_clan_counts = pfam_clan_counts
        self.common_pfam_clans = common_pfam_clans
        
class UniprotRecord(object):
    
    def __init__(self, raw_biopython_record, all_pfam_data):
        
        self.id = raw_biopython_record.id
        self.seq = raw_biopython_record.seq
        self.raw_biopython_record = raw_biopython_record
        
        if all_pfam_data is not None:
            self._set_pfam_domains(all_pfam_data)
        
    def get_closest_pfam_domain(self, location):
        
        pfam_domains = self.pfam_domains.copy()
        
        if len(pfam_domains) == 0:
            return None
            
        pfam_domains['distance'] = pfam_domains.apply(calc_distance_from_pfam_domain, args = (location,), axis = 1)
        return pfam_domains.loc[pfam_domains['distance'].idxmin()]
        
    def __len__(self):
        return len(self.seq)
        
    def __repr__(self):
        return '<UniprotRecord: %s>' % str(self)
        
    def __str__(self):
        return '%s (%d aa long)' % (self.id, len(self.seq))
        
    def _set_pfam_domains(self, all_pfam_data):
        try:
            self.pfam_domains = all_pfam_data.pfam_domains_by_uniprot_id.get_group(self.id)
        except KeyError:
            self.pfam_domains = pd.DataFrame(columns = all_pfam_data.pfam_domains.columns)

def load_pfam_data(config_setup):
    
    pfam_data = pd.read_csv(config_setup.get_path('PFAM_PROTEOME_CSV_FILE_PATH'), sep = r'\t|> <', skiprows = 2, na_values = ['No_clan'], \
            engine = 'python').rename(columns = lambda name: re.sub(r'[ \-]', '_', re.sub(r'[#\<\>]', '', name)).lower())
    pfam_data = pfam_data.rename(columns = {'seq_id': 'uniprot_id'})
    pfam_data['length'] = pfam_data['alignment_end'] - pfam_data['alignment_start'] + 1
    
    pfam_domains = pfam_data[pfam_data['type'] == 'Domain']
    pfam_domains_by_uniprot_id = pfam_domains.groupby('uniprot_id')

    pfam_clan_counts = pfam_domains['clan'].value_counts()
    common_pfam_clans = sorted(set(pfam_clan_counts[pfam_clan_counts >= MIN_COUNTS_FOR_COMMON_CLAN].index))
    
    log(('Loaded %d pfam records, %d of which are domains. Found %d unique clans within these domains, %d of which are common ' + \
            '(with at least %d occurrences).') % (len(pfam_data), len(pfam_domains), len(pfam_clan_counts), len(common_pfam_clans), \
            MIN_COUNTS_FOR_COMMON_CLAN))
    return PfamData(pfam_data, pfam_domains, pfam_domains_by_uniprot_id, pfam_clan_counts, common_pfam_clans)
        
def load_uniprot_data(config_setup, pfam_data = None):
    
    with gzip.open(config_setup.get_path('UNIPROT_PROTEOME_XML_FILE_PATH'), 'r') as f:
        uniprot_records = {record.id: UniprotRecord(record, pfam_data) for record in SeqIO.parse(f, 'uniprot-xml')}

    log('Parsed %d UniProt records.' % len(uniprot_records))
    return uniprot_records

def calc_distance_from_pfam_domain(pfam_record, position):
    
    start = pfam_record['alignment_start']
    end = pfam_record['alignment_end']
    
    if position < start:
        return start - position
    elif position > end:
        return position - end
    else:
        return 0
