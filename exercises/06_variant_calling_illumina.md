# Galaxy for virologist training Exercise 6: Illumina Variant Calling 101
<div class="tables-start"></div>

|**Title**| Galaxy |
|---------|-------------------------------------------|
|**Training dataset:**|  PRJEB43037 - In August 2020 an outbreak of West Nile Virus affected 71 people with meningoencephalitis in Andalusia and 6 more cases in Extremadura (south-west of Spain), causing a total of eight deaths. The virus belonged to the lineage 1 and was relatively similar to previous outbreaks occurred in the Mediterranean region. Here we present a detailed analysis of the outbreak, including an extensive phylogenetic study. This is one of the outbreak samples.
|**Questions:**| <ul><li>What is variant calling?</li><li>What is a vcf file?</li><li>How can I inspect a variant in a bam file to look for false positives?</li><li>How can I make a consensus genome based on a variant calling process?</li></ul>|
|**Objectives**:|<ul><li>Understand variant calling concept</li><li>Learn how to interpret a vcf file</li><li>Learn how to make a reference consensus genome.</li><li>Learn how to visualize mapping and variant calling results</li></ul>|
|**Estimated time**:| 1h |

<div class="tables-end"></div>

## 1. Description
After mapping when we have a re-sequencing experiment the next step usually comprises the variants calling step. Variant calling software tries to identify variants, positions that differ in our reads compared to a reference genome.
Also we may want to have a consensus genome, which we obtain including the variants we just identified in the published reference genome.
We are going to address this type of analysis in this tutorial.

## 2. Upload data to galaxy

### Training dataset
- Experiment info: PRJEB43037, WGS, Illumina MiSeq, paired-end
- Fastq R1: [ERR5310322_1](https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR531/002/ERR5310322/ERR5310322_1.fastq.gz) - url : `ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR531/002/ERR5310322/ERR5310322_1.fastq.gz`
- Fastq R2: [ERR5310322_2](https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR531/002/ERR5310322/ERR5310322_2.fastq.gz)  url : `ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR531/002/ERR5310322/ERR5310322_2.fastq.gz`
- Reference genome NC_009942.1: [fasta](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/875/385/GCF_000875385.1_ViralProj30293/GCF_000875385.1_ViralProj30293_genomic.fna.gz) -- [gff](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/875/385/GCF_000875385.1_ViralProj30293/GCF_000875385.1_ViralProj30293_genomic.gff.gz)

### Create new history
- Click the `+` icon at the top of the history panel and create a new history with the name `mapping 101 tutorial` as explained [here](https://github.com/BU-ISCIII/galaxy_virologist_training/blob/one_week_4day_format/exercises/01_introduction_to_galaxy.md#2-galaxys-history)


### Upload data
Follow the same instructions [here](https://github.com/BU-ISCIII/galaxy_virologist_training/blob/one_week_4day_format/exercises/03_mapping.md#2-upload-data-to-galaxy)

## 3. Preprocess our reads.
Follow instructions [here]()

## 4. Map our reads against our reference genome.
Follow instructions [here](https://github.com/BU-ISCIII/galaxy_virologist_training/blob/one_week_4day_format/exercises/03_mapping.md#map-reads-using-bowtie2)

## 5. Variant Calling.
### Samtools mpileup
1. Search samtools mpileup in the search toolbox.
2. Due to a problem with the last version of mpileup installed in galaxy we need to downgrade to v.2.1.4.
3. Click on Version box in the grey box 
4. Bam files: Bowtie2 bam file
5. Use reference: Use reference/genome from history. NC_009942.1.
6. Set advanced options: Yes
7. Disable read-pair overlap detection: Yes
8. Disable BAQ (per-Base Alignment Quality), see below: Yes
9. max per-file depth; avoids excessive memory usage: 0
10. Minimum base quality for a base to be considered: 20
11. Click execute and wait.

<p align="center"><img src="images/samtools_mpileup_params1.png" alt="samtools mpileup" width="500"></p>
<p align="center"><img src="images/samtools_mpileup_params2.png" alt="samtools mpileup" width="500"></p>
<p align="center"><img src="images/samtools_mpileup_params3.png" alt="samtools mpileup" width="500"></p>

12. Click the :eye: icon on the history and inspect the mpileup output.

### VarScan
1. Search `VarScan Mpileup` in the search toolbox.
2. Samtools pileup dataset: samtools mpileup output
3. Minimum coverage: 10
4. Minimum supporting reads: 5
5. Minimum base quality: 20
6. Minimum variant allele frequency: 0.8
7. Default p-value threshold for calling variants: 0.05
8. Click execute and wait

<p align="center"><img src="images/varscan_params1.png" alt="varscan" width="500"></p>
<p align="center"><img src="images/varscan_params2.png" alt="varscan" width="500"></p>

9. Click the :eye: icon and inspect the vcf file.

### VCF stats
1. Search `bcftools stats` in the search toolbox.
2. VCF/BCF Data: varscan vcf output.
3. Choose a reference genome: use genome/reference from history. Select NC_009942.1.
4. Click execute and wait.
5. Click the :eye: icon and inspect the stats.

<details>
  <summary>How may variants do we have in our vcf file?</summary>
  </br>
  number of SNPs:	463
</details>

## 6. Consensus genome

### Bcftools consensus
1. Search `bcftools consensus` in the search toolbox.
2. VCF/BCF Data: varscan vcf output.
3. Choose a reference genome: use genome/reference from history. Select NC_009942.1.
4. Click execute and wait.

<p align="center"><img src="images/bcftools_consensus_params.png" alt="varscan" width="500"></p>

> Note: for this example we are not going to mask any position with low coverage, this will be addressed in the exercise 8 where a real example will be performed.


