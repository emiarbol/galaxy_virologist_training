# Galaxy for virologist training Exercise 2: Mapping 101

<div class="tables-start"></div>

|**Title**| Galaxy |
|---------|-------------------------------------------|
|**Training dataset:**|  PRJEB43037 - In August 2020 an outbreak of West Nile Virus affected 71 people with meningoencephalitis in Andalusia and 6 more cases in Extremadura (south-west of Spain), causing a total of eight deaths. The virus belonged to the lineage 1 and was relatively similar to previous outbreaks occurred in the Mediterranean region. Here we present a detailed analysis of the outbreak, including an extensive phylogenetic study. This is one of the outbreak samples.
|**Questions:**| <ul><li>What is mapping?</li><li>What is a BAM file?</li><li>Which metrics are important to check after mapping?</ul>|
|**Objectives**:|<ul><li>Understand mapping concept</li><li>Learn how to interpret mapping metrics</li><li>Learn how to visualize mapping mapping results</li></ul>|
|**Estimated time**:| 40 min |

<div class="tables-end"></div>

## 1. Description
One of the most common experiments using massive sequencing are re-sequencing experiments. This type of experiments sequence already known microorganisms where our goal is to discover variation between a aready assembled and known reference and our reads. Mapping is a mandatory step for this kind of experiments, where we need sort all our short sequences(reads) we have in our fastq file without any genomic context.
After the mapping step we are going to transform our fastq file in a bam file where we are going to have information about where a read came from, meaning we are going to have the coordinates where each read is placed in our reference genome.

## 2. Upload data to galaxy

### Training dataset
- Experiment info: PRJEB43037, WGS, Illumina MiSeq, paired-end
- Fastq R1: [ERR5310322_1](https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR531/002/ERR5310322/ERR5310322_1.fastq.gz) - url : `ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR531/002/ERR5310322/ERR5310322_1.fastq.gz`
- Fastq R2: [ERR5310322_2](https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR531/002/ERR5310322/ERR5310322_2.fastq.gz)  url : `ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR531/002/ERR5310322/ERR5310322_2.fastq.gz`
- Reference genome NC_009942.1: [fasta](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/875/385/GCF_000875385.1_ViralProj30293/GCF_000875385.1_ViralProj30293_genomic.fna.gz) -- [gff](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/875/385/GCF_000875385.1_ViralProj30293/GCF_000875385.1_ViralProj30293_genomic.gff.gz)

### Create new history
- Click the `+` icon at the top of the history panel and create a new history with the name `mapping 101 tutorial` as explained [here](https://github.com/BU-ISCIII/galaxy_virologist_training/blob/one_week_4day_format/exercises/01_introduction_to_galaxy.md#2-galaxys-history)


### Upload data
- Import and rename the read files `ERR5310322_1` and `ERR5310322_2`
    1. Click in upload data.
    2. Click in paste/fetch data
    3. Copy url for fastq R1 (select and Ctrl+C) and paste (Ctrl+V).
    4. Click in Start.
    5. Wait until the job finishes (green in history)
    6. Do the same for fastq R2.
<p align="center"><img src="images/upload_data_mapping.png" alt="Upload data mapping" width="900"></p>

- Rename R1 and R2 files.
    1. Click in the pencil in the history for `ERR5310322_1.fastq.gz`
    2. Change the name to `ERR5310322_1`
    3. Do the same for R2.
    
<p align="center"><img src="images/changename1.png" alt="Change name 1" width="900"></p>    

- Import the reference genome.
    
<p align="center"><img src="images/upload_data_mapping2.png" alt="Upload data mapping 2" width="900"></p>

- Rename the reference genome.
    1. Click the pencil for the reference file in the history.
    2. Change the name to `NC_009942.1`

<p align="center"><img src="images/changename2.png" alt="Change name 2" width="900"></p>    

### Map reads using Bowtie2
1. Search bowtie2 software in the search tools box on the left.

<p align="center"><img src="images/search_bowtie2.png" alt="Search bowtie" width="900"></p>   

2. Set bowtie2 parameters:
    - Is this single or paired library: paired.
    - FASTA/Q file #1 : ERR5310322_1
    - FASTA/Q file #2 : ERR5310322_2
    - Will you select a reference genome from your history or use a built-in index? : Use a genome from the history and build index.
    - Do you want to use presets? : Very sensitive local. This setting will hugly affect the mapping results, depending on the dataset/experiment must be tweaked (read [bowtie2 manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml))
    - Save the bowtie2 mapping statistics to the history

<p align="center"><img src="images/bowtie2params1.png" alt="Search bowtie" width="900"></p>   
<p align="center"><img src="images/bowtie2params2.png" alt="Search bowtie" width="900"></p>   

3. Click execute and wait.

### Visualize bam file
