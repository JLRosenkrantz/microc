# Micro-C Workflow

This Snakemake workflow performs **quality control (QC) and processing of Micro-C sequencing data**, following analysis recommendations from [Dovetail Genomics Micro-C Analysis Guide](https://micro-c.readthedocs.io/en/latest/index.html)

***
## Requirements
This workflow requires:
-   **Snakemake** (for workflow automation)
-   **Conda** (for environment management)

***
## Installation & Setup
If you haven't already installed Conda and Snakemake:
1.  **Install Conda**:\
    Follow instructions at Miniconda installation.

2.  **Install Snakemake**:
```bash
    conda install -c conda-forge -c bioconda snakemake
```

3.  **Activate the Snakemake environment** before running the workflow:
```bash
    conda activate snakemake
```

***
## Input Data Requirements

### 1. Raw Data
- Copy raw sequencing files into the **`raw-data/`** directory.
-   Only **`.fastq.gz`** files are accepted.
-   If your files are not .gz compressed, use `gzip` tool to compress

### 2. Sample Naming Convention
- Files must be named using the following format: `Sample1_R1.fastq.gz`, `Sample1_R2.fastq.gz`
-   `_R1` → Read 1 of paired-end sequencing
-   `_R2` → Read 2 of paired-end sequencing\
- **Incorrect names will cause the workflow to fail.**
- **Manually rename files if needed** before running the workflow.

***
## Configuration Setup
Before running the workflow, update the configuration file: `config/config.yml` and set the correct file paths. See example below:

```
reference_genome: "/home/groups/hoolock2/u0/genomes/ucsc/hg38/indexes/bwa/hg38.fa.gz"
genome_file: "/home/groups/hoolock2/u0/genomes/ucsc/hg38/hg38.genome"
chrsizes: "/home/groups/hoolock2/u0/genomes/ucsc/hg38/hg38.chrom.sizes"
threads: 16
```
- **`reference_genome`** → Path to the BWA index file (`.fa.gz`). If this is the first time running for a specific genome, you **must generate index file** (see instructions below).
- **`genome_file`** → Path to the **`.genome`** file. If this is the first time running for a specific genome, you **must generate genome file** (see instructions below).
-   **`chrsizes`** → Path to the chromosome sizes (`.chrom.sizes`) file.
-   **`threads`** → Number of threads to use for multi-threaded tasks.





### Generate BWA index file
Micro-C data is aligned using **Burrows-Wheeler Aligner (BWA)**.\
Before running the workflow, create a **BWA index** for your reference genome.
Run this command in the directory containing the `.fasta` genome sequence:
```
bwa index hg38.fasta
```
This generates the necessary BWA index files.


### Generating a Genome File
A genome file (`.genome`) is required for downstream analysis. If you don't have one, create it from the `.fai` index file:

```
cut -f1,2 hg38.fa.fai > hg38.genome
```

This file contains chromosome names and sizes in **tab-separated format**.

***
## Running the Workflow
Once everything is set up, execute the Snakemake workflow:

### 1. Dry-Run to Check for Issues
Before running, test for missing files or errors:
```
snakemake --use-conda -np
```

### **2️⃣ Run the Workflow**

Execute the full workflow with the desired number of CPU cores:
`snakemake --use-conda --cores 8`


### **3️⃣ Run with Cluster (SLURM)**

If using an HPC with SLURM, submit a job:
`snakemake --use-conda --cores 8 --cluster "sbatch --time=24:00:00 --mem=32G --cpus-per-task=8"`

























**🔹 Running the Workflow**
---------------------------

Once everything is set up, execute the Snakemake workflow:

### **1️⃣ Dry-Run to Check for Issues**

Before running, test for missing files or errors:

bash

CopyEdit

`snakemake --use-conda --cores 1 -n`

### **2️⃣ Run the Workflow**

Execute the full workflow with the desired number of CPU cores:

bash

CopyEdit

`snakemake --use-conda --cores 8`

### **3️⃣ Run with Cluster (SLURM)**

If using an HPC with SLURM, submit a job:

bash

CopyEdit

`snakemake --use-conda --cores 8 --cluster "sbatch --time=24:00:00 --mem=32G --cpus-per-task=8"`

* * * * *

**🔹 Preparing the Genome Index**
---------------------------------

### **1️⃣ BWA Indexing**

Micro-C data is aligned using **Burrows-Wheeler Aligner (BWA)**.\
Before running the workflow, create a **BWA index** for your reference genome.

📌 **Run this command in the directory containing the `.fasta` genome sequence**:

bash

CopyEdit

`bwa index hg38.fasta`

This generates the necessary BWA index files.

### **2️⃣ Generating a Genome File**

A genome file (`.genome`) is required for downstream analysis. If you don't have one, create it from the `.fai` index file:

bash

CopyEdit

`cut -f1,2 /home/groups/hoolock2/u0/genomes/ucsc/hg38/hg38.fa.fai > hg38.genome`

This file contains chromosome names and sizes in **tab-separated format**.

* * * * *

**🔹 Troubleshooting**
----------------------

### 🔸 **Common Errors & Fixes**

| **Issue** | **Solution** |
| --- | --- |
| `Workflow fails due to missing files` | Ensure all paths in `config.yml` are correctly set. |
| `FASTQ file naming error` | Rename files to match the required format (`Sample_R1.fastq.gz`). |
| `BWA index file missing` | Generate it using `bwa index hg38.fasta`. |
| `Out of memory (SLURM jobs fail)` | Increase memory allocation in SLURM (`--mem=64G`). |

* * * * *

**🔹 Citation & Acknowledgments**
---------------------------------

This workflow is adapted from the **Dovetail Micro-C analysis recommendations** and integrates best practices from the **4DN project**.

For additional details, refer to:\
📄 **[Dovetail Micro-C Analysis Guide](https://micro-c.readthedocs.io/en/latest/index.html)**





















# Micro-C Workflow
Snakemake workflow for QC and processing of Micro-C data. Following analysis recommendations from Dovetail (https://micro-c.readthedocs.io/en/latest/index.html)

This workflow requires the use of conda environments and snakemake.
If required download conda and install snakemake.
activate snakemake conda environment before beginning. 

raw files must be copied into raw-data directory
only fastq.gz files are acceptable. Make sure files are .gz compressed.
If not compress using gzip tool.

The sample names are automatically determined based off the fastq.gz file names.  

the files MUST be named as follows Sample1_R1.fastq.gz, Sample1_R2.fastq.gz
It is important that the _R1 and _R2 suffix are used to define read 1 and read 2 of paired end reads. The workflow will fail if named differently.
Manually rename if required.

Before running the snakemake workflow update the config/config.yml file
in this file update correct file paths for the following:
reference_genome: "/home/groups/hoolock2/u0/genomes/ucsc/hg38/indexes/bwa/hg38.fa.gz" this should be the filepath and suffix of the bwa index. you may need to generate this file if this is the first time running for specific genome. See notes below on Genome Index.
genome_file: this should be the path for .genome file - you may need to generate this file if this is the first time running for specific genome. See notes below on Pre-Alignment.
chrsizes: path to .chrom.sizes file
threads: number of threads to run during workflow







Genome Index
In line with the 4DN project guidelines and from our own experience optimal alignment results are obtained with Burrows-Wheeler Aligner (bwa). Prior to alignment, generate a bwa index file for the chosen reference. Run this in directory containing .fasta genome sequence.
example:
bwa index hg38.fasta

Pre-Alignment
For downstream steps you will need a genome file, genome file is a tab delimited file with chromosome names and their respective sizes. If you don’t already have a genome file follow these steps: 
Use the index file to generate the genome file by printing the first two columns into a new file.
example:
cut -f1,2 /home/groups/hoolock2/u0/genomes/ucsc/hg38/hg38.fa.fai > hg38.genome


snakemake --use-conda



Alignment
Now that you have a genome file, index file and a reference fasta file you are all set to align your Micro-C library to the reference. Please note the specific settings that are needed to map mates independently and for optimal results with our proximity library reads.

bwa mem -5SP -T0 -t<threads> <ref.fasta> <MicroC_R1.fastq> <MicroC_R2.fastq> -o <aligned.sam>
Example (one pair of fastq files):

bwa mem -5SP -T0 -t16 hg38.fasta MicroC_2M_R1.fastq MicroC_2M_R2.fastq -o aligned.sam
Example (multiple pairs of fastq files):

bwa mem -5SP -T0 -t16 hg38.fasta <(zcat file1.R1.fastq.gz file2.R1.fastq.gz file3.R1.fastq.gz)



bwa mem -5SP -T0 -t16 /home/groups/hoolock2/u0/genomes/ucsc/hg38/hg38.fa test-data/MicroC_2M_R1.fastq test-data/MicroC_2M_R2.fastq| \
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path /home/groups/hoolock2/u0/genomes/ucsc/hg38/hg38.genome | \
pairtools sort --tmpdir=/home/ubuntu/ebs/temp/ --nproc 16| \
pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats stats.txt| \
pairtools split --nproc-in 8 --nproc-out 8 --output-pairs mapped.pairs --output-sam -| \
samtools view -bS -@16 | \
samtools sort -@16 -o mapped.PT.bam;samtools index mapped.PT.bam









Snakemake workflow for QC and processing of HIC data. Results in on .cool file and TAD calls for the data. All commands used to QC, and process hic-data are contianed in the main Snakefile.

**Prep**:

Place raw HIC reads in a new directory data/raw

```
mkdir -p data/raw
ln -s /path/to/some/data/*.fastq.gz data/raw/
```
Make sure files are compressed with gzip (.fastq.gz NOT .fastq). snakemake workflow will error if files are unzipped .fastq

```
/home/groups/hoolock2/u0/jvc/TAD-RO1/GEO_SUB/nomLeu4_HiC_2_R1.fastq.gz
/home/groups/hoolock2/u0/jvc/TAD-RO1/GEO_SUB/nomLeu4_HiC_2_R2.fastq.gz
/home/groups/hoolock2/u0/bd/tmp_data/gibbon_hic/Vok_NLE_HiC_S3_L006_R1_001.fastq.gz
/home/groups/hoolock2/u0/bd/tmp_data/gibbon_hic/Vok_NLE_HiC_S3_L006_R2_001.fastq.gz

SAMPLES:
    - "Gibbon_nomLeu4"
    - "Gibbon_Vok_NLE"

```


**Configure the file**: src/config.yml

Generate Arima HIC genome restriction site file using the [hicup](https://www.bioinformatics.babraham.ac.uk/projects/hicup/) command [hicup_digester](https://www.bioinformatics.babraham.ac.uk/projects/hicup/) with the flag `--arima` for compatability with the Arima HIC protocol.

All commands used to QC, and process hic-data are contianed in the main Snakefile. The pipeline uses conda for dependency management. Make sure you have installed a recent version of snakemake and conda.
```
example:
hicup_digester --re1 ^GATC,MboI --genome Mouse_mm10 --outdir hicup-digest/ *.fa &
```

**Execution**:

```
snakemake --use-conda -j20
snakemake --use-conda -j20 > 240905h_snakemake.out 2>&1
snakemake --use-conda -j20 > $(date +"%y%m%d%H%M%S")_snakemake.out 2>&1
snakemake --use-conda -j20 cooler_cload zoomify hicFindTADs > $(date +"%y%m%d%H%M%S")_snakemake.out 2>&1


```

**Runtime**

The hicup pipeline is the most resource intensive step that can be expected to run for at least 24 hours for a sample with a sequencing depth of 500 million reads and 8 threads.



##conda activate hicexplorer
hicPlotMatrix --matrix your_data.mcool::resolutions/10000 --log1p --outFileName contact_map.png --title "Gibbon Hi-C Contact Map"

hicPlotMatrix --matrix cool/Gibbon_nomLeu4.mcool::resolutions/10000 --log1p --outFileName Gibbon_nomLeu4_contact_map.png --title "Gibbon_nomLeu4 Hi-C Contact Map"





<!-- 
##Convert .pairs to .hic Using Juicer Tools
cooler dump -t pixels Gibbon_nomLeu4.cool | awk '{print "chr"$1, $2, "chr"$3, $4, $5}' > Gibbon_nomLeu4.pairs &
cooler dump -t bins Gibbon_nomLeu4.cool | cut -f1,3 | sort -u > chrom.sizes &
java -Xmx16g -jar juicer_tools_1.22.01.jar pre Gibbon_nomLeu4.pairs nomLeu4.chrom.sizes Gibbon_nomLeu4.hic &

java -Xmx16g -jar juicer_tools_1.22.01.jar pre -c nomLeu4.chrom.sizes Gibbon_nomLeu4.pairs Gibbon_nomLeu4.hic nomLeu4

pre [options] <infile> <outfile> <genomeID>



java -Xmx16g -jar juicer_tools.3.0.0.jar pre ../pairix/Gibbon_nomLeu4.bsorted.pairs.gz ~/u0/genomes/other/nomLeu4/nomLeu4.chrom.sizes Gibbon_nomLeu4.hic &

mydata.pairs chrom.sizes mydata.hic



hicConvertFormat --matrices Gibbon_nomLeu4.cool --inputFormat cool --output Gibbon_nomLeu4.hic --outputFormat hic &
 -->

<!-- 
hicPlotMatrix --matrix Gibbon_nomLeu4.cool --outFileName HiC_contact_map_Gibbon_nomLeu4.png --log1p --dpi 300 --title "nomLeu4 Genome-wide Hi-C Contact Map"

hicAggregateContacts --matrix Gibbon_nomLeu4.cool --outFileName Gibbon_downsampled.cool --binSize 1000000

hicPlotMatrix --matrix Gibbon_downsampled.cool --outFileName HiC_contact_map_Gibbon_nomLeu4.png --log1p --dpi 300 --title "nomLeu4 Genome-wide Hi-C Contact Map"

cooler zoomify -o Gibbon_nomLeu4.mcool Gibbon_nomLeu4.cool

hicPlotMatrix --matrix Gibbon_nomLeu4.mcool::resolutions/1280000 --outFileName HiC_contact_map_Gibbon_nomLeu4.png --log1p --dpi 300

hicPlotMatrix --matrix Gibbon_nomLeu4.mcool::resolutions/1280000 --outFileName HiC_contact_map_Gibbon_nomLeu4.png --dpi 300

hicPlotMatrix --matrix Gibbon_nomLeu4.cool --outFileName Gibbon_nomLeu4_HiC_contact_map_chr1a.png --log1p --dpi 300 --region chr1a --title "Gibbon_nomLeu4 Hi-C Map: chr1a"

hicPlotMatrix --matrix Gibbon_nomLeu4.cool --outFileName Gibbon_nomLeu4_HiC_contact_map_chr1a.png --dpi 300 --region chr1a --title "Gibbon_nomLeu4 Hi-C Map: chr1a"
 -->















