What is geneffect?
--------

This library combines genomic and proteomic data from various databases (e.g. GENCODE, UCSC's reference genome, UniProt and pfam) into unified gene objects (currently only protein-coding genes). It allows you to infer the functional effects of genetic variations at the gene/protein-level.

If you use geneffect in a work contributing to a scientific publication, we ask that you cite our publication: 

Nadav Brandes, Nathan Linial, Michal Linial, Quantifying gene selection in cancer through protein functional alteration bias, Nucleic Acids Research, gkz546, https://doi.org/10.1093/nar/gkz546


Basic usage
--------

With geneffect installed, you can obtain genomic and proteomic data of protein-coding genes:
    >>> import geneffect
    >>> geneffect_setup = geneffect.Setup('GRCh38')
    [This may take a few minutes...]
    >>> gene_HOXD4, = [gene for gene in geneffect_setup.genes if gene.symbol == 'HOXD4']
    >>> print(gene_HOXD4)
    <Gene: HOXD4, P09016 / <CDSIsoform: ENST00000306324.3 (chr2 (+), 2 CDS exons)>>
    >>> print('The gene %s (%s) is on chromosome %s.' % (gene_HOXD4.symbol, gene_HOXD4.name, gene_HOXD4.canonical_cds_isoform.chromosome))
    The gene HOXD4 (homeobox D4) is on chromosome 2.
    >>> print('DNA sequence: %s' % gene_HOXD4.canonical_cds_isoform.dna_seq)
    DNA sequence: ATGGTCATGAGTTCGTATATGGTGAACTCCAAGTATGTGGACCCCAAGTTCCCTCCGTGCGAGGAGTATTTGCAGGGCGGCTACCTAGGCGAGCAGGGCGCCGACTACTACGGCGGCGGCGCGCAGGGCGCAGACTTCCAGCCCCCGGGGCTCTACCCACGGCCCGACTTCGGTGAGCAGCCTTTCGGAGGCAGCGGCCCCGGGCCTGGCTCGGCGCTGCCTGCGCGGGGTCACGGACAAGAGCCAGGCGGCCCCGGCGGTCACTACGCCGCTCCAGGAGAGCCTTGCCCAGCTCCCCCGGCGCCTCCGCCGGCGCCCCTGCCTGGCGCCCGGGCCTACAGTCAGTCCGACCCCAAGCAGCCGCCCTCCGGGACGGCACTCAAGCAGCCGGCCGTGGTCTACCCCTGGATGAAGAAGGTGCACGTGAATTCGGTGAACCCCAACTACACCGGTGGGGAACCCAAGCGGTCCCGAACGGCCTACACCCGGCAGCAAGTCCTAGAACTGGAAAAAGAATTTCATTTTAACAGGTATCTGACAAGGCGCCGTCGGATTGAAATCGCTCACACCCTGTGTCTGTCGGAGCGCCAGATCAAGATCTGGTTCCAGAACCGGAGGATGAAGTGGAAAAAAGATCATAAGCTGCCCAACACTAAAGGCAGGTCATCGTCCTCATCTTCCTCCTCATCTTGCTCCTCCTCAGTCGCCCCCAGCCAGCATTTACAGCCGATGGCCAAAGACCACCACACGGACCTGACGACCTTA
    >>> print('Protein sequence: %s' % gene_HOXD4.uniprot_record.seq)
    Protein sequence: MVMSSYMVNSKYVDPKFPPCEEYLQGGYLGEQGADYYGGGAQGADFQPPGLYPRPDFGEQPFGGSGPGPGSALPARGHGQEPGGPGGHYAAPGEPCPAPPAPPPAPLPGARAYSQSDPKQPPSGTALKQPAVVYPWMKKVHVNSVNPNYTGGEPKRSRTAYTRQQVLELEKEFHFNRYLTRRRRIEIAHTLCLSERQIKIWFQNRRMKWKKDHKLPNTKGRSSSSSSSSSCSSSVAPSQHLQPMAKDHHTDLTTL
    >>> print('CDS exon coordinates: %s' % ', '.join(['%d..%d' % (exon.chromosome_start, exon.chromosome_end) for exon in gene_HOXD4.canonical_cds_isoform.cds_exons]))
    CDS exon coordinates: 176151634..176152066, 176152608..176152939
    >>> print('UniProt features: %s' % gene_HOXD4.uniprot_record.raw_biopython_record.features)
    UniProt features: [SeqFeature(FeatureLocation(ExactPosition(0), ExactPosition(255)), type='chain', id='PRO_0000200210'), SeqFeature(FeatureLocation(ExactPosition(153), ExactPosition(213)), type='DNA-binding region'), SeqFeature(FeatureLocation(ExactPosition(132), ExactPosition(138)), type='short sequence motif'), SeqFeature(FeatureLocation(ExactPosition(221), ExactPosition(234)), type='compositionally biased region'), SeqFeature(FeatureLocation(ExactPosition(122), ExactPosition(123)), type='sequence variant', id='VAR_067445'), SeqFeature(FeatureLocation(ExactPosition(141), ExactPosition(142)), type='sequence conflict')]

You can also interpret SNPs and their effects on protein-coding genes:
    >>> snp = geneffect_setup.variant_interpreter.process_snp('17', 43082434, 'G', 'C')
    >>> print(snp)
    chr17:43082434G>C [P38398:R1443G]
    >>> print('%s>%s at chr%s:%d' % (snp.ref_dna_seq, snp.alt_dna_seq, snp.chromosome, snp.chromosome_start_coordinate))
    G>C at chr17:43082434
    >>> print('%d CDS gene effects: %s' % (len(snp.cds_gene_effects), snp.cds_gene_effects))
    1 gene effects: [P38398:R1443G]
    >>> snp_cds_gene_effect, = snp.cds_gene_effects
    >>> print(snp_cds_gene_effect)
    P38398:R1443G
    >>> print(snp_cds_gene_effect.affected_gene)
    <Gene: BRCA1, P38398 / <CDSIsoform: ENST00000357654.7 (chr17 (-), 22 CDS exons)>>
    >>> print(snp_cds_gene_effect.is_synonymous(), snp_cds_gene_effect.is_missense(), snp_cds_gene_effect.is_nonsense())
    False True False
    >>> print(snp_cds_gene_effect.protein_coordinate, snp_cds_gene_effect.cds_coordinate, snp_cds_gene_effect.phase)
    1443 4327 0
    >>> print(snp_cds_gene_effect.ref_aa, snp_cds_gene_effect.alt_aa, snp_cds_gene_effect.ref_codon, snp_cds_gene_effect.alt_codon)
    R G CGA GGA
    
As of version 1.2.0, geneffect can also interpret more complex variants:
    >>> variant = geneffect_setup.variant_interpreter.process_variant('17', 41276079, 'A', 'GC') # Frameshift indel
    >>> print(variant)
    chr17:41276079A>GC [KRTAP9-7.382:A->GC [frameshift]]
    >>> print(variant.is_snp(), variant.is_insert(), variant.is_deletion(), variant.is_complex_indel())
    False False False True
    >>> cds_gene_effect, = variant.cds_gene_effects
    >>> print(cds_gene_effect)
    KRTAP9-7.382:A->GC [frameshift]
    >>> print(cds_gene_effect.is_frameshift, cds_gene_effect.phase_change)
    True 1
    >>> variant = geneffect_setup.variant_interpreter.process_variant('17', 41276079, 'A', 'TATAGTG') # Complex amino-acid change
    >>> print(variant)
    chr17:41276079A>TATAGTG [A8MTY7:N128YSD]
    >>> cds_gene_effect, = variant.cds_gene_effects
    >>> print(cds_gene_effect.introduced_stop_codon, cds_gene_effect.destroys_start_codon())
    False False
    >>> print(cds_gene_effect.protein_start_coordinate, cds_gene_effect.ref_protein_seq, cds_gene_effect.alt_protein_seq)
    128 N YSD
    >>> variant = geneffect_setup.variant_interpreter.process_variant('17', 41347060, 'C', 'T') # Canonical splice-site effect
    >>> print(variant)
    chr17:41347060C>TKRT33A: canonical splicing affected in intron #4
    >>> splicing_gene_effect, = variant.splicing_gene_effects
    >>> print(splicing_gene_effect.affected_gene)
    <Gene: KRT33A, O76009 / <CDSIsoform: ENST00000007735.3 (chr17 (-), 7 CDS exons)>>
    >>> print(splicing_gene_effect.affected_intron_index)
    4


Installation
--------

Dependencies:

* numpy
* pandas
* biopython
* interval_tree (https://github.com/moonso/interval_tree)


To install, just run:

    python setup.py install
    
Or:
    
    pip install geneffect
    
(with the latter you will have to copy your configuration file manually)


After installation, you will have to setup your configuration file (by default it is the file .geneffect_config.py in your homedir, or you can define it to be any other file by setting the environment variable GENEFFECT_CONFIG_FILE). The default settings are also available in the file default_config.py within this module.
Just open your configuration file with your favorite editor and follow the instructions within it. In order for this package to work, you will have to download files from five different databases (reference genome from UCSC, gene annotations from GENCODE, metadata of genes from genenames, protein records from UniProt, and, optionally, domain annotations from pfam). The paths of these downloaded files should be set correctly in your configuration file.  
