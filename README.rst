What is geneffect?
--------

This library combines genomic and proteomic data from various databases (e.g. GENCODE, UCSC's reference genome, UniProt and pfam) into unified gene objects (currently only protein-coding genes). It allows you to infer the functional effects of genetic variations at the gene/protein-level.

If you use geneffect in work contributing to a scientific publication, we ask that you cite our publication: 
Nadav Brandes, Nathan Linial, Michal Linial. Modeling Functional Genetic Alteration in Cancer Reveals New Candidate Driver Genes. https://doi.org/10.1101/242354


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
    >>> print('%s>%s at chr%s:%d' % (snp.ref_nt, snp.alt_nt, snp.chromosome, snp.chromosome_coordinate))
    G>C at chr17:43082434
    >>> print('%d gene effects: %s' % (len(snp.gene_effects), snp.gene_effects))
    1 gene effects: [P38398:R1443G]
    >>> snp_gene_effect, = snp.gene_effects
    >>> print(snp_gene_effect)
    P38398:R1443G
    >>> print(snp_gene_effect.affected_gene)
    <Gene: BRCA1, P38398 / <CDSIsoform: ENST00000357654.7 (chr17 (-), 22 CDS exons)>>
    >>> print(snp_gene_effect.is_synonymous(), snp_gene_effect.is_missense(), snp_gene_effect.is_nonsense())
    False True False
    >>> print(snp_gene_effect.protein_coordinate, snp_gene_effect.cds_coordinate, snp_gene_effect.phase)
    1443 4327 0
    >>> print(snp_gene_effect.ref_aa, snp_gene_effect.alt_aa, snp_gene_effect.ref_codon, snp_gene_effect.alt_codon)
    R G CGA GGA


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
