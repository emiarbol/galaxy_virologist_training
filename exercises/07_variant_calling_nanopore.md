# Galaxy for virologist training Exercise 7: Nanopore variant calling and consensus 101

<div class="tables-start"></div>

|**Title**| Galaxy |
|---------|-------------------------------------------|
|**Training dataset:**|  Nanopore Sequencing of a SARS-Cov-2 
|**Questions:**| <ul><li>How can I get a consensus sequence using nanopore reads?</li></ul>|
|**Objectives**:|<ul><li>Use artic protocol</li><li>Learn how to interpret artic results</li></ul>|
|**Estimated time**:| 1 h |

<div class="tables-end"></div>

## 1. Description
Artic network ....blablba

## 2. Upload data to galaxy

### Training dataset
- Experiment info: [sequencing summary](https://github.com/nf-core/test-datasets/blob/viralrecon/nanopore/minion/sequencing_summary.txt)
- fastq: 
    - [fastq1](https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/nanopore/minion/fastq_pass/barcode01/FAO93606_pass_barcode01_7650855b_0.fastq)
    - [fastq2](https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/nanopore/minion/fastq_pass/barcode01/FAO93606_pass_barcode01_7650855b_1.fastq)
    - [fastq3](https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/nanopore/minion/fastq_pass/barcode01/FAO93606_pass_barcode01_7650855b_2.fastq)
- Reference genome MN908947.3 : [fasta](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/858/895/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.fna.gz) --- [gff](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/858/895/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.gff.gz)

### Create new history
- Click the `+` icon at the top of the history panel and create a new history with the name `nanopore assembly 101 tutorial` as explained [here](https://github.com/BU-ISCIII/galaxy_virologist_training/blob/one_week_4day_format/exercises/01_introduction_to_galaxy.md#2-galaxys-history)

### Upload data
- Import and rename the read files `fastq1`, `fastq2`and `fastq3`
> Note: Nanopore reads are commonly splitted in several files that we need to concatenate prior further analysis depending on the software we are going to use.
    1. Click in upload data.
    2. Click in paste/fetch data
    3. Copy url for fastq R1 (select and Ctrl+C) and paste (Ctrl+V).
    4. Click in Start.
    5. Wait until the job finishes (green in history)
    6. Do the same for the remaining files.
<p align="center"><img src="images/upload_data_assemblyNanopore.png" alt="Upload data mapping" width="700"></p>

- Rename files.
    1. Click in the ✏️ in the history for all the files
    2. Change the name to `fastq_X`
    
- Import the reference genome.
    
- Rename the reference genome.
    1. Click the ✏️ for the reference file in the history.
    2. Change the name to `MN908947.3`

### Concatenate reads.
1. Search `Concatenate datasets` using the search toolbox.
2. Select all three fastq files keeping **Cntrl key** clicked.
3. Click execute and wait.
<p align="center"><img src="images/concatenate_params.png" alt="Concatenate" width="700"></p>

## 2. Artic protocol

