geneffect
--------

This library allows you to do the following:
    >>> import geneffect
    >>> geneffect_setup = geneffect.Setup('GRCh38')
    [This may take a few minutes to setup...]
    >>> gene_HOXD4, = [gene for gene in geneffect_setup.genes if gene.symbol == 'HOXD4']
    >>> print(gene_HOXD4)
    <Gene: HOXD4, P09016 / <CDSIsoform: ENST00000306324.3 (chr2 (+), 2 CDS exons)>>
