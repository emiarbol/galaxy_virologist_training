# Galaxy for virologist training
In this training course you will find theory and practice material for introducing yourself to viral genome analysis using Galaxy.

The material includes slides with theory concepts and a bunch of practical exercises, focusing on the interpretation of results.

## Data

- West Nile:
    - Description: In August 2020 an outbreak of West Nile Virus affected 71 people with meningoencephalitis in Andalusia and 6 more cases in Extremadura (south-west of Spain), causing a total of eight deaths. The virus belonged to the lineage 1 and was relatively similar to previous outbreaks occurred in the Mediterranean region. Here we present a detailed analysis of the outbreak, including an extensive phylogenetic study. This is one of the outbreak samples. [paper](https://pubmed.ncbi.nlm.nih.gov/34063166/)
    - Experiment info: PRJEB43037, WGS, Illumina MiSeq, paired-end
    - Fastq R1: [ERR5310322_1](https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR531/002/ERR5310322/ERR5310322_1.fastq.gz)
    - Fastq R2: [ERR5310322_2](https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR531/002/ERR5310322/ERR5310322_2.fastq.gz)
    - Reference genome NC_009942.1: [fasta](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/875/385/GCF_000875385.1_ViralProj30293/GCF_000875385.1_ViralProj30293_genomic.fna.gz) -- [gff](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/875/385/GCF_000875385.1_ViralProj30293/GCF_000875385.1_ViralProj30293_genomic.gff.gz)

- SARS-Cov2 Nanopore
    - Description:Nanopore Sequencing of a SARS-Cov-2 
    - Experiment info: [sequencing summary](https://github.com/nf-core/test-datasets/blob/viralrecon/nanopore/minion/sequencing_summary.txt)
    - fastq: 
        - [fastq1](https://github.com/nf-core/test-datasets/blob/viralrecon/nanopore/minion/fastq_pass/barcode01/FAO93606_pass_barcode01_7650855b_0.fastq)
        - [fastq2](https://github.com/nf-core/test-datasets/blob/viralrecon/nanopore/minion/fastq_pass/barcode01/FAO93606_pass_barcode01_7650855b_1.fastq)
        - [fastq3](https://github.com/nf-core/test-datasets/blob/viralrecon/nanopore/minion/fastq_pass/barcode01/FAO93606_pass_barcode01_7650855b_2.fastq)
    - Reference genome MN908947.3 : [fasta](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/858/895/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.fna.gz) --- [gff](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/858/895/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.gff.gz)

- SARS-Cov2 Illumina:
    - Description:
    - Experiment info:
    - Fastq R1:
    - Fastq R2:
    - Reference genome MN908947.3 : [fasta](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/858/895/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.fna.gz) --- [gff](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/858/895/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.gff.gz)

## Slides
### Day 1

- **Talk 1:** [Uso de la secuenciación masiva en virología](slides/XXXXXX.pdf)
- **Talk 2:** [Introducción teórica a Galaxy](slides/02_introducción_galaxy.pdf)
- [**Exercise 1**](exercises/01_introduction_to_galaxy.md) -- [Download pdf](exercises/01_introduction_to_galaxy.pdf)

### Day 2

- **Talk 3:** [Análisis de calidad y preprocesamiento de los datos genómicos](slides/Day2/03_quality_preprocessing.pdf)
- **Talk 4:** [Fase de mapping](slides/curso_ViralGalaxy_session_Mapping.pdf)
- **Talk 5:** [Ensamblado de genomas](slides/XX.pdf)
- [**Exercise 2**](exercises/02_quality.md) -- [Download pdf](exercises/02_quality.pdf)
- [**Exercise 3**](exercises/03_mapping.md) -- [Download pdf](exercises/03_mapping.pdf)
- [**Exercise 4**](exercises/04_assembly_illumina.md) -- [Download pdf](exercises/04_assembly_illumina.pdf)
- [**Exercise 5**](exercises/05_assembly_nanopore.md) -- [Download pdf](exercises/05_assembly_nanopore.pdf)


### Day 3

- **Talk 6:** [Llamada a variantes: identificación de SNPs y obtención de secuencia consenso](slides/curso_ViralGalaxy_session_VariantCallingConsensus.pdf)
- [**Exercise 6**](exercises/06_variant_calling_illumina.md) -- [Download pdf](exercises/06_variant_calling_illumina.pdf)

### Day 4

- **Talk 7:** [Viralrecon](slides/viralrecon.pdf)
- [**Exercise 8**](exercises/04_viralrecon.md) -- [Download pdf](exercises/04_viralrecon.pdf)
