# Micro-C Workflow
Snakemake workflow for QC and processing of Micro-C data. Following analysis recommendations from Dovetail (https://micro-c.readthedocs.io/en/latest/index.html)

activate snakemake conda environment before beginning

raw files must be copied into raw-data directory
only fastq.gz is acceptable

the files MUST be named as follows sampleID_readX.fastq.gz
for example Sample1_R1.fastq.gz, Sample1_R2.fastq.gz

snakemake --use-conda


reference in a fasta file format, e.g. hg38

REFFASTA="/home/groups/hoolock2/u0/genomes/ucsc/hg38/hg38.fa"

Pre-Alignment
For downstream steps you will need a genome file, genome file is a tab delimited file with chromosome names and their respective sizes. If you don’t already have a genome file follow these steps: 
Use the index file to generate the genome file by printing the first two columns into a new file.
example:
cut -f1,2 /home/groups/hoolock2/u0/genomes/ucsc/hg38/hg38.fa.fai > hg38.genome


In line with the 4DN project guidelines and from our own experience optimal alignment results are obtained with Burrows-Wheeler Aligner (bwa). Prior to alignment, generate a bwa index file for the chosen reference. Run this in directory containing .fasta genome sequence.
example:
bwa index hg38.fasta


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















