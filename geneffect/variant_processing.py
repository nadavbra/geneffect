from __future__ import absolute_import, division, print_function

from collections import defaultdict

# Install from: https://github.com/moonso/interval_tree
from interval_tree import IntervalTree

from .util import as_biopython_seq, codon_to_aa

class SNP(object):

    def __init__(self, chromosome, chromosome_coordinate, ref_nt, alt_nt, gene_effects):
        self.chromosome = chromosome
        self.chromosome_coordinate = chromosome_coordinate
        self.ref_nt = ref_nt
        self.alt_nt = alt_nt
        self.gene_effects = gene_effects
        
    def __repr__(self):
        
        effects_repr = ', '.join(map(str, self.gene_effects))
        
        if len(self.gene_effects) > 0:
            effects_repr = ' [%s]' % effects_repr
        
        return 'chr%s:%d%s>%s%s' % (self.chromosome, self.chromosome_coordinate, self.ref_nt, self.alt_nt, effects_repr)
        
class SNPGeneEffect(object):
    
    def __init__(self, affected_gene, affected_cds_exon, cds_coordinate, protein_coordinate, ref_codon, alt_codon, phase, \
            ref_aa, alt_aa):
        
        self.affected_gene = affected_gene
        self.affected_cds_exon = affected_cds_exon
        self.cds_coordinate = cds_coordinate
        self.protein_coordinate = protein_coordinate
        self.ref_codon = ref_codon
        self.alt_codon = alt_codon
        self.phase = phase
        self.ref_aa = ref_aa
        self.alt_aa = alt_aa
        
    def is_nonsense(self):
        return self.alt_aa == '*'
        
    def is_synonymous(self):
        return self.ref_aa == self.alt_aa
        
    def is_missense(self):
        return not self.is_nonsense() and not self.is_synonymous()
        
    def __repr__(self):
        return '%s:%s%d%s' % (self.affected_gene.uniprot_record.id, self.ref_aa, self.protein_coordinate, self.alt_aa)

class VariantInterpreter(object):
    
    def __init__(self, genome_reader, genes):
        
        self.genome_reader = genome_reader
        genes_per_chromosome = defaultdict(list)

        for gene in genes:
            genes_per_chromosome[gene.canonical_cds_isoform.chromosome].append(gene)
            
        self.chromosome_interpreters = {chromosome: ChromosomeVariantInterpreter(chr_name, chr_genes, genome_reader) for chr_name, chr_genes in \
                genes_per_chromosome.items()}
                
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
        
    def process_snp(self, coordinate, ref_nt, alt_nt):
        
        ref_nt = as_biopython_seq(ref_nt)
        alt_nt = as_biopython_seq(alt_nt)
        
        if self.chromosome in self.genome_reader:
            read_ref_nt = self.genome_reader.read_seq(self.chromosome, coordinate, coordinate)
            assert ref_nt == read_ref_nt, 'Read %s instead of %s at chr%s:%d' % (read_ref_nt, ref_nt, self.chromosome, coordinate)
        
        gene_effects = []
        
        for gene in self.get_genes_in_locus(coordinate):
            for cds_exon in gene.canonical_cds_isoform.cds_exons:
                if cds_exon.chromosome_start <= coordinate <= cds_exon.chromosome_end:
                    
                    if gene.canonical_cds_isoform.strand == '+':
                        cds_ref_nt = ref_nt
                        cds_alt_nt = alt_nt
                        cds_coordinate = cds_exon.isoform_start + coordinate - cds_exon.chromosome_start
                    else:
                        cds_ref_nt = ref_nt.reverse_complement()
                        cds_alt_nt = alt_nt.reverse_complement()
                        cds_coordinate = cds_exon.isoform_start + cds_exon.chromosome_end - coordinate
                    
                    gene_effect = build_snp_gene_effect(gene, cds_exon, cds_coordinate, cds_ref_nt, cds_alt_nt)
                    
                    if gene_effect is not None:
                        gene_effects.append(gene_effect)
                    
                    break
                    
        return SNP(self.chromosome, coordinate, ref_nt, alt_nt, gene_effects)
    
    def get_genes_in_locus(self, coordinate):
        if self.gene_interval_tree is None:
            return []
        else:
            return self.gene_interval_tree.find_range([coordinate, coordinate])
        
def build_snp_gene_effect(gene, cds_exon, cds_coordinate, cds_ref_nt, cds_alt_nt):

    protein_coordinate = (cds_coordinate - 1) // 3 + 1
    phase = (cds_coordinate - 1) % 3
    
    ref_codon = gene.canonical_cds_isoform.dna_seq[(3 * (protein_coordinate - 1)):(3 * protein_coordinate)]
    
    # It turns out that in version GRCh38 of the reference genome, we end up with some genes having a canonical_cds_isoform
    # of length which is not a multiply of 3. Therefore, nt substitution in the end of these genes result a partial codon.
    # While I'm not sure why this happens (even though these genes do end up giving the correct aa sequence), this small patch
    # should avoid the problem.
    if len(ref_codon) != 3:
        return None
    
    assert ref_codon[phase] == cds_ref_nt
    alt_codon = ref_codon[:phase] + cds_alt_nt + ref_codon[(phase + 1):]

    ref_aa = codon_to_aa(ref_codon, gene.canonical_cds_isoform.codon_table)
    assert gene.uniprot_record.seq[protein_coordinate - 1] == ref_aa
    alt_aa = codon_to_aa(alt_codon, gene.canonical_cds_isoform.codon_table)
    
    return SNPGeneEffect(gene, cds_exon, cds_coordinate, protein_coordinate, ref_codon, alt_codon, phase, ref_aa, alt_aa)
        
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
