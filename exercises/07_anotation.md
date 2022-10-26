# Galaxy for virologist training Exercise 7: Illumina Variant Annotation 101
<div class="tables-start"></div>

|**Title**| Galaxy |
|---------|-------------------------------------------|
|**Training dataset:**|  PRJEB43037 - In August 2020, an outbreak of West Nile Virus affected 71 people with meningoencephalitis in Andalusia and 6 more cases in Extremadura (south-west of Spain), causing a total of eight deaths. The virus belonged to the lineage 1 and was relatively similar to previous outbreaks occurred in the Mediterranean region. Here, we present a detailed analysis of the outbreak, including an extensive phylogenetic study. This is one of the outbreak samples.
|**Questions:**| <ul><li>Which effects have variants in the genome?</li></ul>|
|**Objectives**:|<ul><li>Understand the importance of variants effect significance.</li></ul>|
|**Estimated time**:| 1h |

<div class="tables-end"></div>

## 1. Description
After performing variant calling, we want to know which is the importance of the variants in the viral genome. In order to give sense to the variants, we need to know in which gene they are, and which are their effects.

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
Follow instructions [here](https://github.com/BU-ISCIII/galaxy_virologist_training/edit/one_week_format/exercises/06_variant_calling_illumina.md#5-variant-calling)

## 6. Variants annotation
### Load annotation file for West Nile genome.
1. Load gff file
2. Upload file
3. Paste/Fetch data: [gff](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/875/385/GCF_000875385.1_ViralProj30293/GCF_000875385.1_ViralProj30293_genomic.gff.gz)
4. Rename dataset: NC_009942.1.

### Snpeff build
1. Search `snpeff build` in the search toolbox.
2. Name of the database: WestNile.
3. Input annotations are in: GFF
4. GFF dataset to build database from: NC_009942.1 gff
5. Choose the source of the reference genome: History. NC_009942.1 fasta.
6. Click execute and wait.

<p align="center"><img src="images/snpeff_build_params1.png" alt="snpeff build" width="500"></p>

### Snpeff eff
1. Search `snpeff eff` in the search toolbox.
2. Sequence changes (SNPs, MNPs, InDels): ivar vcf file
3. Genome source: Custom snpEff database in your history. Snpeff build output.
4. Create CSV report, useful for downstream analysis (-csvStats): Yes.
5. Click execute and wait.

<p align="center"><img src="images/snpeff_eff_params1.png" alt="snpeff eff" width="500"></p>

6. Click the :eye: icon in the SnpEff html output and check the results.

### SnpSift: transfrom vcf snpeff to table.
1.  Search `SnpSift ExtractFields` in the search toolbox.
2.  Variant input file in VCF format: snpeff eff vcf output.
3.  Fields to extract: CHROM POS ID REF ALT FILTER ANN[*].EFFECT ANN[*].GENE ANN[*].FEATURE ANN[*].HGVS_C ANN[*].HGVS_P
4.  One effect per line: Yes.
5.  Click execute and wait.
6.  Click the :eye: icon in the snpsift output and check the results.

<p align="center"><img src="images/snpsift_params1.png" alt="snpsift" width="500"></p>
