
# Bash Script for ITS Sequencing Data Analysis

This script performs several steps in the analysis of ITS sequencing data.

Please note that this script requires certain tools like `fastqc`, `multiqc`, and `qiime`. Make sure you have them installed and properly configured in your environment.

## Step 1: Download the ITS sequencing data from SRA database
This step uses the `fastq-dump` command from the SRA Toolkit to download the sequencing data from the SRA database. The `--split-files` option is used to split paired-end reads into separate files, and the `--gzip` option is used to compress the output files using gzip.

```bash
fastq-dump --split-files --gzip
ERR6377054	ERR6377055	ERR6377056
ERR6377058	ERR6377059	ERR6377060
ERR6377062	ERR6377137	ERR6377138
ERR6377139	ERR6377141	ERR6377142
ERR6377143	ERR6377144	ERR6377145
```

## Step 2: Quality Check

This step involves checking the quality of the downloaded sequencing data. The `fastqc` command is used to generate quality control reports for each of the fastq files. The `multiqc` command is then used to aggregate these reports into a single report.

```bash
fastqc *.gz
multiqc .
```

## Step 3: Generate a manifest.txt and formatting for qiime2
This step involves creating a manifest file that lists all the fastq files along with their absolute paths. This file will be used in the next step for importing the data into QIIME 2.

```bash
ls -d "$PWD"/* > manifest.txt
```

## Step 4: Pack the paired-end data
This step involves importing the paired-end sequencing data into QIIME 2. The `qiime tools import` command is used for this purpose. The `--type 'SampleData[PairedEndSequencesWithQuality]'` option specifies the type of data being imported. The `--input-path manifest.csv` option specifies the path to the manifest file created in the previous step. The `--output-path demux` option specifies the output directory where the imported data will be stored. The `--input-format PairedEndFastqManifestPhred33` option specifies the format of the input data.

```bash
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.csv \
  --output-path demux \
  --input-format PairedEndFastqManifestPhred33
```

## Step 5: Sequences quality visualisation (before the trimming)
This step involves visualizing the quality of the imported sequences before trimming. The `qiime demux summarize` command is used for this purpose. The `--i-data demux.qza` option specifies the input data, and the `--o-visualization demux` option specifies the output visualization.

```bash
qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux
```

## Step 6: Cutadapt for the adaptors and primers
This step involves trimming adaptors and primers from sequences using Cutadapt, a tool designed to clean biological sequences, especially high-throughput sequencing reads. This removes unwanted biases in your data before downstream analysis.

```bash
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux.qza \
  --p-front-f TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGAACCWGCGGARGGATCA \
  --p-front-r GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGCTGCGTTCTTCATCGATGC \
  --p-error-rate 0 \
  --o-trimmed-sequences trimmed-demux.qza \
  --verbose > primer_trimming.log
```

## Step 7: Sequences quality visualisation (after the trimming)
This step involves visualizing sequence quality after trimming adaptors and primers. This allows you to check if trimming was successful and if further preprocessing steps are necessary.

```bash
qiime demux summarize \
  --i-data trimmed-demux.qza \
  --o-visualization trimmed-demux
```

## Step 8: Denoise Paired Sequences
This step involves denoising paired-end sequences using DADA2, a pipeline for detecting and correcting (or removing) Illumina amplicon sequence data. This includes model-based quality filtering, dereplication (removal of duplicates), sample inference, merging of paired-end reads, and chimera removal using DADA2. Due to the high degree of variation in length across different species,   `--p-trim-left-f 0`  and `--p-trim-left-r 0` need to be used to improve the biological reads.

```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimmed-demux_its.qza \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-max-ee-f 4 \
  --p-max-ee-r 4 \
  --p-n-threads 2 \
  --o-representative-sequences rep-seqs-its.qza \
  --o-table table-its.qza \
  --o-denoising-stats stats-its.qza
```

## Step 9: Tabulate Sequences
This step involves tabulating sequence data for easy visualization and inspection. It provides an overview of sequence variants in your sample.

```bash
qiime feature-table tabulate-seqs \
  --i-data rep-seqs-its.qza \
  --o-visualization rep-seqs-its.qzv
```

## Step 10: Classify Sequences with UNITE Classifier
This step involves classifying sequences using a pre-trained Naive Bayes classifier and the UNITE database, which is specifically designed for fungal ITS sequences.

```bash
qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/shared/scratch/pyau/qiime2-ref/unite-9-dynamic-s-all-29.11.2022-Q2-2023.5.qza \
  --i-reads rep-seqs-its.qza \
  --o-classification taxonomy-its.qza
```

## Step 11: Tabulate Metadata
This step involves tabulating metadata for easy visualization. It provides an overview of taxonomic assignments in your sample.

```bash
qiime metadata tabulate \
  --m-input-file taxonomy-its.qza \
  --o-visualization taxonomy-its.qzv
```

## Step 12: Taxonomy-based Filtering of Tables and Sequences
This step involves filtering tables and sequences based on taxonomy. Sequences that are classified as mitochondria or chloroplasts are excluded from further analysis as they do not represent microbial diversity.

```bash
qiime taxa filter-table \
--i-table table-its.qza \
--i-taxonomy taxonomy-its.qza \
--p-include p__ \
--p-exclude mitochondria,chloroplast \
--o-filtered-table table-its-with-phyla-no-mitochondria-no-chloroplast.qza
  
qiime taxa filter-seqs \
--i-sequences rep-seqs-its.qza \
--i-taxonomy taxonomy-its.qza \
--p-include p__ \
--p-exclude mitochondria,chloroplast \
--o-filtered-sequences rep-seqs-its-with-phyla-no-mitochondria-no-chloroplast.qza
```
