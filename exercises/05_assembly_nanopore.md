# Galaxy for virologist training Exercise 4: Illumina Assembly 101

<div class="tables-start"></div>

|**Title**| Galaxy |
|---------|-------------------------------------------|
|**Training dataset:**|  PRJEB43037 - In August 2020 an outbreak of West Nile Virus affected 71 people with meningoencephalitis in Andalusia and 6 more cases in Extremadura (south-west of Spain), causing a total of eight deaths. The virus belonged to the lineage 1 and was relatively similar to previous outbreaks occurred in the Mediterranean region. Here we present a detailed analysis of the outbreak, including an extensive phylogenetic study. This is one of the outbreak samples.
|**Questions:**| <ul><li>What is assembly?</li><li>How can I evaluate my assembly?</li></ul>|
|**Objectives**:|<ul><li>Understand assembly concept</li><li>Learn how to interpret assembly quality control metrics</li></ul>|
|**Estimated time**:| 40 min |

<div class="tables-end"></div>

## 1. Description
When we don't have a reference genome to map against it, or when we don't want to have any bias in the genome reconstruction what we need to do is a de novo assembly. This type of analysis tries to reconstruct the original genome without any template using only the reads. 
Some considerations:
- When we assembly as longer the reads and as longer the size of the library fragments the easier it gets for the assembler. That's why pacbio or nanopore are recommended for assembly.
- It's almost imposible to reconstruct the entire genome of a large microorganism with only one sequencing, although it can be done for small viruses.
- Assembly is not recommended for amplicon based libraries due to the depth of coverage uneveness and the amplicons intrinsic bias.

## 2. Upload data to galaxy

### Training dataset
- Experiment info: PRJEB43037, WGS, Illumina MiSeq, paired-end
- Fastq R1: [ERR5310322_1](https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR531/002/ERR5310322/ERR5310322_1.fastq.gz) - url : `ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR531/002/ERR5310322/ERR5310322_1.fastq.gz`
- Fastq R2: [ERR5310322_2](https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR531/002/ERR5310322/ERR5310322_2.fastq.gz)  url : `ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR531/002/ERR5310322/ERR5310322_2.fastq.gz`
- Reference genome NC_009942.1: [fasta](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/875/385/GCF_000875385.1_ViralProj30293/GCF_000875385.1_ViralProj30293_genomic.fna.gz) -- [gff](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/875/385/GCF_000875385.1_ViralProj30293/GCF_000875385.1_ViralProj30293_genomic.gff.gz)

### Create new history
- Click the `+` icon at the top of the history panel and create a new history with the name `illumina assembly 101 tutorial` as explained [here](https://github.com/BU-ISCIII/galaxy_virologist_training/blob/one_week_4day_format/exercises/01_introduction_to_galaxy.md#2-galaxys-history)


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
    1. Click in the ✏️ in the history for `ERR5310322_1.fastq.gz`
    2. Change the name to `ERR5310322_1`
    3. Do the same for R2.
    
<p align="center"><img src="images/changename1.png" alt="Change name 1" width="900"></p>    

- Import the reference genome.
    
<p align="center"><img src="images/upload_data_mapping2.png" alt="Upload data mapping 2" width="900"></p>

- Rename the reference genome.
    1. Click the ✏️ for the reference file in the history.
    2. Change the name to `NC_009942.1`

<p align="center"><img src="images/changename2.png" alt="Change name 2" width="900"></p>    

### Assemble reads with Spades
