from __future__ import absolute_import, division, print_function

from operator import itemgetter
import io
import json

import numpy as np
import pandas as pd

from .util import log, append_if_not_none, get_unique_or_none, as_biopython_seq

class Gene(object):
    
    def __init__(self, symbol, name, refseq_ids, uniprot_record, canonical_cds_isoform):
        self.symbol = symbol
        self.name = name
        self.refseq_ids = refseq_ids
        self.uniprot_record = uniprot_record
        self.canonical_cds_isoform = canonical_cds_isoform
        
    def get_identifier(self):
    
        identifiers = []
    
        if self.symbol is not None:
            identifiers.append(self.symbol)
            
        identifiers.append(self.uniprot_record.id)
        return ', '.join(identifiers)
        
    def __repr__(self):
        return '<Gene: %s / %s>' % (self.get_identifier(), self.canonical_cds_isoform) 

class CDSIsoform(object):
    
    def __init__(self, transcript_id, chromosome, strand, cds_chromosomal_coordinates, genome_reader):
        self.transcript_id = transcript_id
        self.chromosome = chromosome
        self.strand = strand
        self.codon_table = 'Vertebrate Mitochondrial' if self.chromosome == 'M' else 'Standard'
        self._set_cds_exons(cds_chromosomal_coordinates, genome_reader)
        self.dna_seq = as_biopython_seq(''.join([str(cds_exon.seq) for cds_exon in self.cds_exons]))
        self.translated_seq = self.dna_seq.translate(table = self.codon_table)
        
    def __len__(self):
        return len(self.dna_seq)
        
    def __repr__(self):
        return '<CDSIsoform: %s (chr%s (%s), %d CDS exons)>' % (self.transcript_id, self.chromosome, self.strand, \
                len(self.cds_exons))
                
    def _set_cds_exons(self, cds_chromosomal_coordinates, genome_reader):
        
        self.cds_exons = []
        last_isoform_end = 0
        
        for chromosome_start, chromosome_end in cds_chromosomal_coordinates:
            cds_exon = CDSExon(self.chromosome, self.strand, chromosome_start, chromosome_end, last_isoform_end + 1, genome_reader)
            self.cds_exons.append(cds_exon)
            last_isoform_end = cds_exon.isoform_end
        
class CDSExon(object):

    def __init__(self, chromosome, strand, chromosome_start, chromosome_end, isoform_start, genome_reader):
        
        self.chromosome = chromosome
        self.strand = strand
        
        self.chromosome_start = chromosome_start
        self.chromosome_end = chromosome_end
        self.length = chromosome_end - chromosome_start + 1
        
        self.isoform_start = isoform_start
        self.isoform_end = isoform_start + self.length - 1
        self.phase = (self.isoform_start - 1) % 3
        
        self.positive_strand_seq = genome_reader.read_seq(chromosome, chromosome_start, chromosome_end)
        self.seq = self.positive_strand_seq if strand == '+' else self.positive_strand_seq.reverse_complement()
        
    def __repr__(self):
        return 'CDS Exon %d-%d at chr%s:%d-%d (%s)' % (self.isoform_start, self.isoform_end, self.chromosome, self.chromosome_start, \
                self.chromosome_end, self.strand)

def load_genes(config_setup, genome_reader, uniprot_records):

    genes = []
    cds_records = _load_cds_records(config_setup)

    for uniprot_id, uniprot_cds_rows in cds_records.groupby('uniprot_id'):
        if uniprot_id in uniprot_records:
            append_if_not_none(genes, _parse_gene(uniprot_records[uniprot_id], uniprot_cds_rows, genome_reader))
            
    log('Parsed %d genes.' % len(genes))
    return genes
                
def _parse_gene(uniprot_record, cds_rows, genome_reader):
        
    cds_isoforms = [_parse_cds_isoform(transcript_id, transcript_cds_rows, genome_reader) for transcript_id, transcript_cds_rows in \
            cds_rows.groupby('transcript_id')]
    matching_cds_isoforms = [isoform for isoform in cds_isoforms if isoform.translated_seq == uniprot_record.seq]
            
    if len(matching_cds_isoforms) > 0:
        symbol = get_unique_or_none(cds_rows['symbol'].unique())
        name = get_unique_or_none(cds_rows['name'].unique())
        refseq_ids = set.union(*map(set, cds_rows['refseqs']))
        return Gene(symbol, name, refseq_ids, uniprot_record, matching_cds_isoforms[0])
    else:
        return None     

def _parse_cds_isoform(transcript_id, transcript_cds_rows, genome_reader):
    
    transcript_cds_rows = transcript_cds_rows.sort_values('exon_number')
    exon_numbers = transcript_cds_rows['exon_number']
    assert list(exon_numbers) == list(range(exon_numbers.min(), exon_numbers.max() + 1))
    
    chromosome, = transcript_cds_rows['chr'].unique()
    strand, = transcript_cds_rows['strand'].unique()
    cds_coordinates = [(row['start'], row['end']) for _, row in transcript_cds_rows.iterrows()]
    return CDSIsoform(transcript_id, chromosome, strand, cds_coordinates, genome_reader)
    
def _load_cds_records(config_setup):
    
    gene_annotations = _load_gene_annotations(config_setup)
    gene_meta_data = _load_gene_meta_data(config_setup)
    
    cds_records = gene_annotations[gene_annotations['type'] == 'CDS']
    cds_records = pd.merge(cds_records, gene_meta_data, left_on = 'gene_id', right_on = 'ensembel_id')

    cds_records['exon_number'] = cds_records['extra_fields'].apply(itemgetter('exon_number')).astype(int)
    cds_records['chr'] = cds_records['chr'].apply(lambda name: name[3:])
    cds_records['transcript_id'] = cds_records['extra_fields'].apply(itemgetter('transcript_id'))
    
    log('Constructed %d CDS records.' % len(cds_records))
    return cds_records
    
def _load_gene_annotations(config_setup):
    
    gene_annotations = pd.read_csv(config_setup.get_path('GENCODE_GENE_ANNOTATIONS_CSV_FILE_PATH'), sep = '\t', comment = '#', \
            names = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'extra_fields'])
    gene_annotations['extra_fields'] = gene_annotations['extra_fields'].apply(_parse_annotation_extra_fields)
    gene_annotations['gene_id'] = gene_annotations['extra_fields'].apply(lambda extra_fields: extra_fields['gene_id'].split('.')[0])
    gene_annotations['gene_type'] = gene_annotations['extra_fields'].apply(itemgetter('gene_type'))
    
    log('Loaded %d gene annotations.' % len(gene_annotations))
    return gene_annotations
    
def _load_gene_meta_data(config_setup):
    
    with io.open(config_setup.get_path('GENENAMES_GENE_META_DATA_JSON_FILE_PATH'), 'r', encoding = 'utf-8') as f:
        raw_meta_data = json.load(f)
    
    ensembel_ids = []
    uniprot_ids = []
    refseqs = []
    symbols = []
    names = []
    gene_families = []

    for record in raw_meta_data['response']['docs']:
        for uniprot_id in record.get('uniprot_ids', []):
            ensembel_ids += [record.get('ensembl_gene_id', np.nan)]
            uniprot_ids += [uniprot_id]
            refseqs += [record.get('refseq_accession', [])]
            symbols += [record['symbol']]
            names += [record['name']]
            gene_families += [record.get('gene_family', [])]
            
    gene_meta_data = pd.DataFrame({'ensembel_id': ensembel_ids, 'uniprot_id': uniprot_ids, 'refseqs': refseqs, \
            'symbol': symbols, 'name': names, 'gene_families': gene_families})
    log('Loaded %d gene meta data records.' % len(gene_meta_data))
    return gene_meta_data
    
def _parse_annotation_extra_fields(raw_extra_fields):

    extra_fields = {}

    for raw_extra_field in raw_extra_fields[:-1].split(';'):
        key, raw_value = raw_extra_field.strip().split(' ')
        value = raw_value.strip('"')
        extra_fields[key] = value
        
    return extra_fields
