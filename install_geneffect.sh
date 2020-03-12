#! /bin/tcsh -ef


# Determine the Python version to use.

# If a the environment variable PYTHON_CMD exists, then just use it.
if ($?PYTHON_CMD) then
    set python_cmd = $PYTHON_CMD
    goto determine_final_python_version
endif

ask_default_python_version:
printf "Do you want to install geneffect in your default Python version (`python --version |& cat`)? Yes/No (Default: Yes): "
set response = $<
set response = `echo $response | tr "[:upper:]" "[:lower:]"`

if ($response == "" || $response == "y" || $response == "yes") then
    set python_cmd = "python"
else if ($response == "n" || $response == "no") then
    printf "Please specify an alternative version of Python: "
    set python_cmd = $<
else
    echo "Please specify either Yes or No."
    goto ask_default_python_version
endif

determine_final_python_version:
echo "Will install geneffect for the Python version `$python_cmd --version |& cat` ($python_cmd)"


# Determine a temp working directory.

# If a the environment variable PYTHON_CMD exists, then just use it.
if ($?TEMP_DIR) then
    set temp_dir = $TEMP_DIR
    goto determine_final_temp_dir
endif

ask_default_temp_dir:
printf "Do you want to use /tmp as a temporary working directory? Yes/No (Default: Yes): "
set response = $<
set response = `echo $response | tr "[:upper:]" "[:lower:]"`

if ($response == "" || $response == "y" || $response == "yes") then
    set temp_dir = "/tmp"
else if ($response == "n" || $response == "no") then
    printf "Please specify an alternative directory path: "
    set temp_dir = $<
else
    echo "Please specify either Yes or No."
    goto ask_default_temp_dir
endif

determine_final_temp_dir:
set temp_dir = `echo $temp_dir | sed 's:/*$::'`
echo "Will use $temp_dir as a temporary working directory."


# Determine the data directory.

# If a the environment variable DATA_DIR exists, then just use it.
if ($?DATA_DIR) then
    set data_dir = $DATA_DIR
    goto determine_final_data_dir
endif

ask_data_dir:
printf "Do you want to use ~/data as the directory for all the data files required by geneffect? Yes/No (Default: Yes): "
set response = $<
set response = `echo $response | tr "[:upper:]" "[:lower:]"`

if ($response == "" || $response == "y" || $response == "yes") then
    set data_dir = "~/data"
else if ($response == "n" || $response == "no") then
    printf "Please specify an alternative directory path: "
    set data_dir = $<
else
    echo "Please specify either Yes or No."
    goto ask_data_dir
endif

determine_final_data_dir:
set data_dir = `echo $data_dir | sed 's:/*$::'`
echo "Will use $data_dir as the directory for data files."


# Install dependencies.

echo "Installing dependencies..."
$python_cmd -m pip install numpy pandas biopython interval_tree 


# Cloning and installing geneffect.

mkdir -p ${temp_dir}/geneffect
git clone https://github.com/nadavbra/geneffect.git ${temp_dir}/geneffect
cd ${temp_dir}/geneffect
$python_cmd ./setup.py install
cp ./geneffect/default_config.py ~/.geneffect_config.py
cd -
rm -fr ${temp_dir}/geneffect


# Downloading the data files and update the configuration.

sed -i "/DATA_DIR = /c\DATA_DIR = os.path.expanduser('${data_dir}')" ~/.geneffect_config.py

echo "Downloading the hg19 human reference genome..."
mkdir -p ${data_dir}/human_reference_genome/hg19
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz -O ${data_dir}/human_reference_genome/hg19/chr1.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr2.fa.gz -O ${data_dir}/human_reference_genome/hg19/chr2.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr3.fa.gz -O ${data_dir}/human_reference_genome/hg19/chr3.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr4.fa.gz -O ${data_dir}/human_reference_genome/hg19/chr4.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr5.fa.gz -O ${data_dir}/human_reference_genome/hg19/chr5.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr6.fa.gz -O ${data_dir}/human_reference_genome/hg19/chr6.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr7.fa.gz -O ${data_dir}/human_reference_genome/hg19/chr7.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr8.fa.gz -O ${data_dir}/human_reference_genome/hg19/chr8.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr9.fa.gz -O ${data_dir}/human_reference_genome/hg19/chr9.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr10.fa.gz -O ${data_dir}/human_reference_genome/hg19/chr10.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr11.fa.gz -O ${data_dir}/human_reference_genome/hg19/chr11.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr12.fa.gz -O ${data_dir}/human_reference_genome/hg19/chr12.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr13.fa.gz -O ${data_dir}/human_reference_genome/hg19/chr13.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr14.fa.gz -O ${data_dir}/human_reference_genome/hg19/chr14.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr15.fa.gz -O ${data_dir}/human_reference_genome/hg19/chr15.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr16.fa.gz -O ${data_dir}/human_reference_genome/hg19/chr16.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr17.fa.gz -O ${data_dir}/human_reference_genome/hg19/chr17.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr18.fa.gz -O ${data_dir}/human_reference_genome/hg19/chr18.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr19.fa.gz -O ${data_dir}/human_reference_genome/hg19/chr19.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr20.fa.gz -O ${data_dir}/human_reference_genome/hg19/chr20.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr21.fa.gz -O ${data_dir}/human_reference_genome/hg19/chr21.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz -O ${data_dir}/human_reference_genome/hg19/chr22.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chrX.fa.gz -O ${data_dir}/human_reference_genome/hg19/chrX.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chrY.fa.gz -O ${data_dir}/human_reference_genome/hg19/chrY.fa.gz
set DOLLAR = '$'
wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log${DOLLAR}=seqview&db=nuccore&report=fasta&sort=&id=251831106&from=begin&to=end&maxplex=1" -O ${data_dir}/human_reference_genome/hg19/chrM.fa
gunzip ${data_dir}/human_reference_genome/hg19/*.gz

echo "Downloading the GRCh38 human reference genome..."
mkdir -p ${data_dir}/human_reference_genome/hg38
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz -O ${data_dir}/human_reference_genome/hg38/chr1.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr2.fa.gz -O ${data_dir}/human_reference_genome/hg38/chr2.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr3.fa.gz -O ${data_dir}/human_reference_genome/hg38/chr3.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr4.fa.gz -O ${data_dir}/human_reference_genome/hg38/chr4.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr5.fa.gz -O ${data_dir}/human_reference_genome/hg38/chr5.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr6.fa.gz -O ${data_dir}/human_reference_genome/hg38/chr6.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr7.fa.gz -O ${data_dir}/human_reference_genome/hg38/chr7.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr8.fa.gz -O ${data_dir}/human_reference_genome/hg38/chr8.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr9.fa.gz -O ${data_dir}/human_reference_genome/hg38/chr9.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr10.fa.gz -O ${data_dir}/human_reference_genome/hg38/chr10.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr11.fa.gz -O ${data_dir}/human_reference_genome/hg38/chr11.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr12.fa.gz -O ${data_dir}/human_reference_genome/hg38/chr12.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr13.fa.gz -O ${data_dir}/human_reference_genome/hg38/chr13.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr14.fa.gz -O ${data_dir}/human_reference_genome/hg38/chr14.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr15.fa.gz -O ${data_dir}/human_reference_genome/hg38/chr15.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr16.fa.gz -O ${data_dir}/human_reference_genome/hg38/chr16.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr17.fa.gz -O ${data_dir}/human_reference_genome/hg38/chr17.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr18.fa.gz -O ${data_dir}/human_reference_genome/hg38/chr18.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr19.fa.gz -O ${data_dir}/human_reference_genome/hg38/chr19.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz -O ${data_dir}/human_reference_genome/hg38/chr20.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz -O ${data_dir}/human_reference_genome/hg38/chr21.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz -O ${data_dir}/human_reference_genome/hg38/chr22.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chrX.fa.gz -O ${data_dir}/human_reference_genome/hg38/chrX.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chrY.fa.gz -O ${data_dir}/human_reference_genome/hg38/chrY.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chrM.fa.gz -O ${data_dir}/human_reference_genome/hg38/chrM.fa.gz
gunzip ${data_dir}/human_reference_genome/hg38/*.gz

echo "Downloading the GENCODE files..."
mkdir -p ${data_dir}/gencode
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz -O ${data_dir}/gencode/gencode.v19.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz -O ${data_dir}/gencode/gencode.v26.annotation.gtf.gz

echo "Downloading the meta data of the human genome from genenames..."
mkdir -p ${data_dir}/genenames
wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/json/non_alt_loci_set.json -O ${data_dir}/genenames/non_alt_loci_set.json

echo "Downloading the human proteome from UniProt..."
mkdir -p ${data_dir}/uniprot
wget "http://www.uniprot.org/uniprot/?sort=score&desc=&compress=yes&query=organism:%22Homo%20sapiens%20[9606]%22%20AND%20reviewed:yes&fil=&format=xml&force=yes" -O ${data_dir}/uniprot/uniprot_human_reviewed.xml.gz

echo "Downloading all human Pfam domains (relase 30.0)..."
mkdir -p ${data_dir}/pfam
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/proteomes/9606.tsv.gz -O ${data_dir}/pfam/9606.tsv.gz

echo "Succefully installed geneffect."
