'''
Instructions:

Two versions of reference genomes are currently configured: hg19 and GRCh38 (both for human). In order to support more reference genome
versions (potentially of other organisms as well, although this was not tested), you first need to add them to the REF_GENOME_NAMES
dictionary (by providing all of the alternative names that the user will be able to use to configure the library to work with those
versions). 

For each supported version of the reference genome, all the variables in this file need to be configured as well:
* REFERENCE_GENOME_DIR (the directory with the reference genome sequences)
* GENCODE_GENE_ANNOTATIONS_CSV_FILE_PATH (the file with GENCODE's gene annotations)
* GENENAMES_GENE_META_DATA_JSON_FILE_PATH (the file with genenames' data)
* UNIPROT_PROTEOME_XML_FILE_PATH (the file with UniProt's entire proteome)
* PFAM_PROTEOME_CSV_FILE_PATH (the file with pfam's entire proteome)

Each of these variables is a list of entries. Each entry is a tuple of two elements. The first defines the relevant versions
of the reference genome; the second defines a dir/file path for the relevant data (downloaded from one of the required resources).
In order to obtain the data files, follow the instructions in the comments within each of these variables.  Feel free to modify
the paths to whatever locations are convenient for you (or to the paths where the required resources are already downloaded in your
system).
'''

from __future__ import absolute_import, division, print_function

import os

DATA_DIR = os.path.join(os.path.expanduser('~'), 'data')

REF_GENOME_NAMES = {
    # hg19 (human)
    'hg19': 'hg19',
    'GRCh37': 'hg19',
    # GRCh38 (human)
    'hg38': 'GRCh38',
    'GRCh38': 'GRCh38',
}

REFERENCE_GENOME_DIR = [
    # The reference genome sequences of all human chromosomes (chrXXX.fa.gz files) can be downloaded from UCSC's FTP site at:
    # ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/ (for version hg19)
    # ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/ (for version hg38/GRCh38)
    # The chrXXX.fa.gz files need to be uncompressed to obtain chrXXX.fa files.
    # IMPORTANT: In version hg19 there's an inconsistency in the reference genome of the M chromosome between UCSC and RegSeq/GENCODE,
    # so the file chrM.fa has to be taken from RefSeq (NC_012920.1) instead of UCSC, from:
    # https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&sort=&id=251831106&from=begin&to=end&maxplex=1
    # In GRCh38 all downloaded files should remain they are.
    (('hg19',), os.path.join(DATA_DIR, 'human_reference_genome/hg19')),
    (('GRCh38',), os.path.join(DATA_DIR, 'human_reference_genome/hg38')),
]

GENCODE_GENE_ANNOTATIONS_CSV_FILE_PATH = [
    # GENCODE's gene annotations for the entire human genome in version hg19 can be downloaded from their FTP site at:
    # ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
    (('hg19',), os.path.join(DATA_DIR, 'gencode/gencode.v19.annotation.gtf.gz')),
    # For version GRCh38:
    # ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz
    (('GRCh38',), os.path.join(DATA_DIR, 'gencode/gencode.v26.annotation.gtf.gz')),
]

GENENAMES_GENE_META_DATA_JSON_FILE_PATH = [
    # Meta data (including different names) of all human genes can be downloaded from genenames' FTP site at:
    # ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/json/non_alt_loci_set.json
    (('hg19', 'GRCh38'), os.path.join(DATA_DIR, 'genenames/non_alt_loci_set.json')),
]

UNIPROT_PROTEOME_XML_FILE_PATH = [
    # All reviewed human proteins can be downloaded in XML format from the following query page:
    # http://www.uniprot.org/uniprot/?sort=score&desc=&compress=yes&query=organism:%22Homo%20sapiens%20[9606]%22%20AND%20reviewed:yes&fil=&format=xml&force=yes
    (('hg19', 'GRCh38'), os.path.join(DATA_DIR, 'uniprot/uniprot_human_reviewed.xml.gz')),
]

PFAM_PROTEOME_CSV_FILE_PATH = [
    # Domain annotations for the entire human proteome can be downloaded from pfam's FTP site at:
    # ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/proteomes/9606.tsv.gz
    (('hg19', 'GRCh38'), os.path.join(DATA_DIR, 'pfam/9606.tsv.gz')),
]
