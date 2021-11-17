# Galaxy for virologist training Exercise 2: Mapping 101

<div class="tables-start"></div>

|**Title**| Galaxy |
|---------|-------------------------------------------|
|**Training dataset:**|  PRJEB43037 - In August 2020 an outbreak of West Nile Virus affected 71 people with meningoencephalitis in Andalusia and 6 more cases in Extremadura (south-west of Spain), causing a total of eight deaths. The virus belonged to the lineage 1 and was relatively similar to previous outbreaks occurred in the Mediterranean region. Here we present a detailed analysis of the outbreak, including an extensive phylogenetic study. This is one of the outbreak samples.
|**Questions:**| <ul><li>What is mapping?</li><li>What is a BAM file?</li><li>Which metrics are important to check after mapping?</ul>|
|**Objectives**:|<ul><li>Understand mapping concept</li><li>Learn how to interpret mapping metrics</li><li>Learn how to visualize mapping mapping results</li></ul>|
|**Estimated time**:| 1h 45 min |

<div class="tables-end"></div>

## 1. Description
One of the most common experiments using massive sequencing are re-sequencing experiments. This type of experiments sequence already known microorganisms where our goal is to discover variation between a aready assembled and known reference and our reads. Mapping is a mandatory step for this kind of experiments, where we need sort all our short sequences(reads) we have in our fastq file without any genomic context.
After the mapping step we are going to transform our fastq file in a bam file where we are going to have information about where a read came from, meaning we are going to have the coordinates where each read is placed in our reference genome.

## 2. Upload data to galaxy

### Training dataset
- Experiment info: PRJEB43037, WGS, Illumina MiSeq, paired-end
- Fastq R1: [ERR5310322_1](https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR531/002/ERR5310322/ERR5310322_1.fastq.gz)
- Fastq R2: [ERR5310322_2](https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR531/002/ERR5310322/ERR5310322_2.fastq.gz)
- Reference genome NC_009942.1: [fasta](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/875/385/GCF_000875385.1_ViralProj30293/GCF_000875385.1_ViralProj30293_genomic.fna.gz) -- [gff](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/875/385/GCF_000875385.1_ViralProj30293/GCF_000875385.1_ViralProj30293_genomic.gff.gz)

### Create new history
- Click the `+` icon at the top of the history panel and create a new history with the name `mapping 101 tutorial` as explained [here](https://github.com/BU-ISCIII/galaxy_virologist_training/blob/one_week_4day_format/exercises/01_introduction_to_galaxy.md#2-galaxys-history)
- Import and rename the read files `ERR5310322_1` and `ERR5310322_2`
<p align="center"><img src="images/upload_data_mapping.png" alt="Upload data mapping" width="900"></p>

- Import the reference genome.
- Rename the reference genome.
