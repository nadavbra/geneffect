from __future__ import absolute_import, division, print_function

from collections import defaultdict

# Install from: https://github.com/moonso/interval_tree
from interval_tree import IntervalTree

from .util import as_biopython_seq, join_seqs

class VariantAtTheEndOfInvalidGeneException(Exception):
    '''
    gene_data_loader ends up with some genes having a canonical_cds_isoform of length which is not a multiply of 3. Therefore, variants at the
    end of these genes result partial codons. I'm not sure why this happens (even though these genes do end up giving the correct aa sequence),
    so this exception is thrown whenever these rare occasions are encountered.
    '''
    pass

class VariantInterpreter(object):
    
    def __init__(self, genome_reader, genes):
        
        self.genome_reader = genome_reader
        genes_per_chromosome = defaultdict(list)

        for gene in genes:
            genes_per_chromosome[gene.canonical_cds_isoform.chromosome].append(gene)
            
        self.chromosome_interpreters = {chr_name: ChromosomeVariantInterpreter(chr_name, chr_genes, genome_reader) for chr_name, chr_genes in \
                genes_per_chromosome.items()}
                
    def process_variant(self, chromosome, start_coordinate, ref_dna_seq, alt_dna_seq):
        return self.get_chromsome_interpreter(chromosome).process_variant(start_coordinate, ref_dna_seq, alt_dna_seq)
        
    def process_snp(self, chromosome, coordinate, ref_nt, alt_nt):
        return self.get_chromsome_interpreter(chromosome).process_snp(coordinate, ref_nt, alt_nt)
        
    def get_chromsome_interpreter(self, chromosome):
    
        if chromosome == 'MT':
            chromosome = 'M'
    
        if chromosome in self.chromosome_interpreters:
            return self.chromosome_interpreters[chromosome]
        else:
            return ChromosomeVariantInterpreter(chromosome, [], self.genome_reader)
        
class ChromosomeVariantInterpreter(object):

    def __init__(self, chromosome, genes, genome_reader):
        self.chromosome = chromosome
        self.genes = genes
        self.genome_reader = genome_reader
        self.gene_interval_tree = _build_gene_interval_trees(genes)
        
    def process_variant(self, start_coordinate, ref_dna_seq, alt_dna_seq):
        
        ref_dna_seq = as_biopython_seq(ref_dna_seq)
        alt_dna_seq = as_biopython_seq(alt_dna_seq)
        end_coordinate = start_coordinate + len(ref_dna_seq) - 1
        
        if self.chromosome in self.genome_reader:
            read_ref_dna_seq = self.genome_reader.read_seq(self.chromosome, start_coordinate, end_coordinate)
            assert ref_dna_seq == read_ref_dna_seq, 'Read %s instead of %s at chr%s:%d' % (read_ref_dna_seq, ref_dna_seq, self.chromosome, \
                    start_coordinate)
        
        cds_gene_effects = []
        splicing_gene_effects = []
        
        for gene in self.get_genes_in_locus(start_coordinate, end_coordinate):
            
            try:
            
                cds_gene_effect = process_cds_gene_effect(gene, start_coordinate, ref_dna_seq, alt_dna_seq)
            
                if cds_gene_effect is not None:
                    cds_gene_effects.append(cds_gene_effect)
            except VariantAtTheEndOfInvalidGeneException:
                pass
                
            splicing_gene_effects.extend(process_splicing_gene_effects(gene, start_coordinate, end_coordinate))
        
        if len(ref_dna_seq) == 1 and len(alt_dna_seq) == 1:
            return SNP(self.chromosome, start_coordinate, ref_dna_seq, alt_dna_seq, cds_gene_effects, splicing_gene_effects)
        else:
            return Variant(self.chromosome, start_coordinate, ref_dna_seq, alt_dna_seq, cds_gene_effects, splicing_gene_effects)
            
    def process_snp(self, coordinate, ref_nt, alt_nt):
        assert len(ref_nt) == 1 and len(alt_nt) == 1
        snp = self.process_variant(coordinate, ref_nt, alt_nt)
        assert isinstance(snp, SNP)
        return snp
            
    def get_genes_in_locus(self, start_coordinate, end_coordinate):
        if self.gene_interval_tree is None:
            return []
        else:
            return self.gene_interval_tree.find_range([start_coordinate, end_coordinate])

class Variant(object):

    def __init__(self, chromosome, chromosome_start_coordinate, ref_dna_seq, alt_dna_seq, cds_gene_effects, splicing_gene_effects):
        
        self.chromosome = chromosome
        self.cds_gene_effects = cds_gene_effects
        self.splicing_gene_effects = splicing_gene_effects
        
        self.ref_dna_seq = ref_dna_seq
        self.alt_dna_seq = alt_dna_seq
        self.ref_dna_len = len(ref_dna_seq)
        self.alt_dna_len = len(alt_dna_seq)
        self.diff_dna_len = self.alt_dna_len - self.ref_dna_len
        assert self.ref_dna_len > 0 or self.alt_dna_len > 0
        
        self.chromosome_start_coordinate = chromosome_start_coordinate
        self.chromosome_end_coordinate = chromosome_start_coordinate + self.ref_dna_len - 1
        
    def is_snp(self):
        return self.ref_dna_len == 1 and self.alt_dna_len == 1
        
    def is_insert(self):
        return self.ref_dna_len == 0
        
    def is_deletion(self):
        return self.alt_dna_len == 0
        
    def is_complex_indel(self):
        return not self.is_snp() and not self.is_insert() and not self.is_deletion()
        
    def get_canonical_splicing_affected_genes(self):
        
        canonical_splicing_affected_genes = set()
        
        for splicing_gene_effect in self.splicing_gene_effects:
            assert isinstance(splicing_gene_effect, CanonicalSplicingGeneEffect)
            canonical_splicing_affected_genes.add(splicing_gene_effect.affected_gene)
        
        return canonical_splicing_affected_genes
        
    def __repr__(self):
        
        effects_repr = ', '.join(map(str, self.cds_gene_effects + self.splicing_gene_effects))
        
        if len(self.cds_gene_effects) > 0:
            effects_repr = ' [%s]' % effects_repr
        
        return 'chr%s:%d%s>%s%s' % (self.chromosome, self.chromosome_start_coordinate, self.ref_dna_seq, self.alt_dna_seq, effects_repr)

class SNP(Variant):
    def __init__(self, chromosome, chromosome_coordinate, ref_nt, alt_nt, snp_cds_gene_effects, splicing_gene_effects):
        assert all([isinstance(snp_cds_gene_effect, SnpCDSGeneEffect) for snp_cds_gene_effect in snp_cds_gene_effects])
        super(SNP, self).__init__(chromosome, chromosome_coordinate, ref_nt, alt_nt, snp_cds_gene_effects, splicing_gene_effects)
        assert self.is_snp()
        
class CanonicalSplicingGeneEffect(object):
    
    def __init__(self, affected_gene, affected_intron_index):
        self.affected_gene = affected_gene
        self.affected_intron_index = affected_intron_index
        
    def __repr__(self):
        return '%s: canonical splicing affected in intron #%d' % (self.affected_gene.symbol, self.affected_intron_index)

class CDSGeneEffect(object):
    
    def __init__(self, affected_gene, affected_cds_exons, cds_start_coordinate, ref_dna_seq, alt_dna_seq):
        
        self.affected_gene = affected_gene
        self.affected_cds_exons = affected_cds_exons
        
        self.ref_dna_seq = ref_dna_seq
        self.alt_dna_seq = alt_dna_seq
        self.ref_dna_len = len(ref_dna_seq)
        self.alt_dna_len = len(alt_dna_seq)
        assert self.ref_dna_len > 0 or self.alt_dna_len > 0
        
        self.cds_start_coordinate = cds_start_coordinate
        self.cds_end_coordinate = cds_start_coordinate + self.ref_dna_len - 1
        assert self.affected_gene.canonical_cds_isoform.dna_seq[(self.cds_start_coordinate - 1):self.cds_end_coordinate] == ref_dna_seq
        self.full_alt_dna_seq = self.affected_gene.canonical_cds_isoform.dna_seq[:(self.cds_start_coordinate - 1)] + self.alt_dna_seq + \
                self.affected_gene.canonical_cds_isoform.dna_seq[(self.cds_end_coordinate - 1):]
        
        self.diff_dna_len = self.alt_dna_len - self.ref_dna_len
        self.phase_change = self.diff_dna_len % 3
        self.is_phased = (self.phase_change == 0)
        self.is_frameshift = not self.is_phased
        
    def is_snp(self):
        return self.ref_dna_len == 1 and self.alt_dna_len == 1
        
    def is_insert(self):
        return self.ref_dna_len == 0
        
    def is_deletion(self):
        return self.alt_dna_len == 0
        
    def is_complex_indel(self):
        return not self.is_snp() and not self.is_insert() and not self.is_deletion()
        
    def destroys_start_codon(self):
        return self.full_alt_dna_seq[:3] != _DNA_START_CODON
        
    def __repr__(self):
        return '%s.%d:%s->%s%s' % (self.affected_gene.symbol, self.cds_start_coordinate, self.ref_dna_seq, self.alt_dna_seq, ' [frameshift]' if \
                self.is_frameshift else '')
        
class PhasedCDSGeneEffect(CDSGeneEffect):

    def __init__(self, affected_gene, affected_cds_exons, cds_start_coordinate, ref_dna_seq, alt_dna_seq):
        
        super(PhasedCDSGeneEffect, self).__init__(affected_gene, affected_cds_exons, cds_start_coordinate, ref_dna_seq, alt_dna_seq)
        assert self.is_phased
        
        if len(affected_gene.canonical_cds_isoform.dna_seq) % 3 != 0 and self.cds_end_coordinate > 3 * (len(affected_gene.canonical_cds_isoform.dna_seq) // 3):
            # The gene's CDS length is not a multiply of 3, and this variant touches the last (partial) codon.
            raise VariantAtTheEndOfInvalidGeneException()
        
        self.start_phase = (self.cds_start_coordinate - 1) % 3
        self.end_phase = (self.cds_end_coordinate - 1) % 3
        
        self.protein_start_coordinate = (self.cds_start_coordinate - 1) // 3 + 1
        self.protein_end_coordinate = (self.cds_end_coordinate - 1) // 3 + 1
        self.ref_protein_len = self.protein_end_coordinate - self.protein_start_coordinate + 1
        self.ref_codons = self.affected_gene.canonical_cds_isoform.dna_seq[(3 * (self.protein_start_coordinate - 1)):(3 * self.protein_end_coordinate)]
        self.ref_protein_seq = self.ref_codons.translate(table = self.affected_gene.canonical_cds_isoform.codon_table)
        assert len(self.ref_protein_seq) == self.ref_protein_len
        assert '*' not in self.ref_protein_seq
        
        self.alt_codons = self.ref_codons[:self.start_phase] + self.alt_dna_seq + ('' if self.end_phase >= 2 else self.ref_codons[(self.end_phase - 2):])
        assert len(self.alt_codons) % 3 == 0
        self.alt_protein_seq = self.alt_codons.translate(table = self.affected_gene.canonical_cds_isoform.codon_table)
        assert len(self.alt_codons) == 3 * len(self.alt_protein_seq)
        self.introduced_stop_codon = ('*' in self.alt_protein_seq)
        
        if self.introduced_stop_codon:
            self.alt_protein_seq = self.alt_protein_seq[:(self.alt_protein_seq.find('*') + 1)]
            self.alt_protein_len = len(self.alt_protein_seq) - 1
            self.diff_protein_len = len(self.affected_gene.canonical_cds_isoform.translated_seq) - self.protein_start_coordinate - self.alt_protein_len + 1
        else:
            self.alt_protein_len = len(self.alt_protein_seq)
            self.diff_protein_len = self.alt_protein_len - self.ref_protein_len
        
    def __repr__(self):
        return '%s:%s%d%s' % (self.affected_gene.uniprot_record.id, self.ref_protein_seq, self.protein_start_coordinate, self.alt_protein_seq)

class SnpCDSGeneEffect(PhasedCDSGeneEffect):
    
    def __init__(self, affected_gene, affected_cds_exon, cds_coordinate, cds_ref_nt, cds_alt_nt):
    
        super(SnpCDSGeneEffect, self).__init__(affected_gene, [affected_cds_exon], cds_coordinate, cds_ref_nt, cds_alt_nt)
        assert self.is_snp()
        
        self.affected_cds_exon = affected_cds_exon
        self.cds_coordinate = cds_coordinate
        self.cds_ref_nt = cds_ref_nt
        self.cds_alt_nt = cds_alt_nt
        
        assert self.cds_start_coordinate == self.cds_end_coordinate
        self.cds_coordinate = self.cds_start_coordinate
        assert self.protein_start_coordinate == self.protein_end_coordinate
        self.protein_coordinate = self.protein_start_coordinate
        
        assert len(self.ref_codons) == 3
        self.ref_codon = self.ref_codons
        assert len(self.alt_codons) == 3
        self.alt_codon = self.alt_codons
        
        assert self.start_phase == self.end_phase
        self.phase = self.start_phase
        
        assert len(self.ref_protein_seq) == 1
        self.ref_aa = self.ref_protein_seq
        assert len(self.alt_protein_seq) == 1
        self.alt_aa = self.alt_protein_seq
        
    def is_nonsense(self):
        return self.alt_aa == '*'
        
    def is_synonymous(self):
        return self.ref_aa == self.alt_aa
        
    def is_missense(self):
        return not self.is_nonsense() and not self.is_synonymous()
        
    def __repr__(self):
        return '%s:%s%d%s' % (self.affected_gene.uniprot_record.id, self.ref_aa, self.protein_coordinate, self.alt_aa)
        
def process_cds_gene_effect(gene, start_chromosome_coordinate, ref_dna_seq, alt_dna_seq):

    end_chromosome_coordinate = start_chromosome_coordinate + len(ref_dna_seq) - 1

    affected_cds_exons = []
    cds_start_coordinate = None
    cds_ref_dna_seqs = []
    cds_alt_dna_seq = None

    for cds_exon in gene.canonical_cds_isoform.cds_exons:
    
        overlapping_start = max(start_chromosome_coordinate, cds_exon.chromosome_start)
        overlapping_end = min(end_chromosome_coordinate, cds_exon.chromosome_end)
        
        if overlapping_start <= overlapping_end:
        
            affected_cds_exons.append(cds_exon)
            cds_overlapping_ref_dna_seq = ref_dna_seq[(overlapping_start - start_chromosome_coordinate):(overlapping_end - start_chromosome_coordinate + 1)]
            
            if gene.canonical_cds_isoform.strand == '+':
                cds_overlapping_start = cds_exon.isoform_start + overlapping_start - cds_exon.chromosome_start
                cds_overlapping_end = cds_exon.isoform_start + overlapping_end - cds_exon.chromosome_start
            else:
                cds_overlapping_ref_dna_seq = cds_overlapping_ref_dna_seq.reverse_complement()
                cds_overlapping_start = cds_exon.isoform_start + cds_exon.chromosome_end - overlapping_end
                cds_overlapping_end = cds_exon.isoform_start + cds_exon.chromosome_end - overlapping_start
                
            assert gene.canonical_cds_isoform.dna_seq[(cds_overlapping_start - 1):cds_overlapping_end] == cds_overlapping_ref_dna_seq
            
            if cds_start_coordinate is None:
                
                cds_start_coordinate = cds_overlapping_start
                
                if gene.canonical_cds_isoform.strand == '+':
                    cds_alt_dna_seq = alt_dna_seq[(overlapping_start - start_chromosome_coordinate):]
                else:
                    cds_alt_dna_seq = alt_dna_seq.reverse_complement()[(overlapping_start - start_chromosome_coordinate):]
            
            cds_ref_dna_seqs.append(cds_overlapping_ref_dna_seq)
        
    if len(affected_cds_exons) == 0:
        return None
    else:
        
        cds_ref_dna_seq = join_seqs(cds_ref_dna_seqs)
        
        if len(cds_ref_dna_seq) == 1 and len(cds_alt_dna_seq) == 1:
            affected_cds_exon, = affected_cds_exons
            return SnpCDSGeneEffect(gene, affected_cds_exon, cds_start_coordinate, cds_ref_dna_seq, cds_alt_dna_seq)
        elif len(cds_ref_dna_seq) % 3 == len(cds_alt_dna_seq) % 3:
            return PhasedCDSGeneEffect(gene, affected_cds_exons, cds_start_coordinate, cds_ref_dna_seq, cds_alt_dna_seq)
        else:
            return CDSGeneEffect(gene, affected_cds_exons, cds_start_coordinate, cds_ref_dna_seq, cds_alt_dna_seq)
            
def process_splicing_gene_effects(gene, start_chromosome_coordinate, end_chromosome_coordinate):

    for intron_index in range(1, len(gene.canonical_cds_isoform.cds_exons)):
    
        exon_before = gene.canonical_cds_isoform.cds_exons[intron_index - 1]
        exon_after = gene.canonical_cds_isoform.cds_exons[intron_index]
    
        if gene.canonical_cds_isoform.strand == '-':
            exon_before, exon_after = exon_after, exon_before
            
        intron_start_chromosome_coordinate = exon_before.chromosome_end + 1
        intron_end_chromosome_coordinate = exon_after.chromosome_start - 1
            
        if _overlaps(start_chromosome_coordinate, end_chromosome_coordinate, intron_start_chromosome_coordinate, intron_start_chromosome_coordinate + 1) or \
                _overlaps(start_chromosome_coordinate, end_chromosome_coordinate, intron_end_chromosome_coordinate - 1, intron_end_chromosome_coordinate):
                
            # The first or last 2 nucleotides of the intron are hit by the variant.
            yield CanonicalSplicingGeneEffect(gene, intron_index)
        
def _build_gene_interval_trees(genes):

    if len(genes) == 0:
        return None

    segments = []
    max_coordinate = 1
    
    for gene in genes:
        start, end = _get_gene_locus(gene)
        segments += [(start, end, gene)]
        max_coordinate = max(max_coordinate, end)
        
    return IntervalTree(segments, 1, max_coordinate)
        
def _get_gene_locus(gene):
    
    all_coordinates = set()
    
    for cds_exon in gene.canonical_cds_isoform.cds_exons:
        all_coordinates.add(cds_exon.chromosome_start)
        all_coordinates.add(cds_exon.chromosome_end)
        
    return min(all_coordinates), max(all_coordinates)
    
def _overlaps(start1, end1, start2, end2):
    return max(start1, start2) <= min(end1, end2)

_DNA_START_CODON = as_biopython_seq('ATG')
