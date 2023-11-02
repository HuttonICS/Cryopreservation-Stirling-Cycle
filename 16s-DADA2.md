
# Bash Script for 16S Sequencing Data Analysis

This script performs several steps in the analysis of 16S sequencing data.

Please note that this script requires certain tools like `fastqc`, `multiqc`, and `qiime`. Make sure you have them installed and properly configured in your environment. Also, please replace the paths with the actual paths where your files are located.

## Step 1: Download the 16S sequencing data from SRA database

This step uses the `fastq-dump` command from the SRA Toolkit to download the sequencing data from the SRA database. The `--split-files` option is used to split paired-end reads into separate files, and the `--gzip` option is used to compress the output files using gzip.

```bash
fastq-dump --split-files --gzip ERR6375874	ERR6375957	ERR6375997	ERR6376032	ERR6376063	ERR6376114	ERR6376145	ERR6376173	ERR6376240	ERR6376291	ERR6376334	ERR6376495	ERR6376576	ERR6376581	ERR6376586	ERR6376591	ERR6376592	ERR6376812	ERR6376813	ERR6376814	ERR6376815	ERR6377043	ERR6377045	ERR6377046	ERR6377048	ERR6377049	ERR6377050	ERR6377051	ERR6377052	ERR6377054	ERR6377055	ERR6377056	ERR6377058	ERR6377059	ERR6377060	ERR6377062	ERR6377137	ERR6377138 ERR6377139 ERR6377141 ERR6377142 ERR6377143 ERR6377144 ERR6377145
```

## Step 2: Quality Check

This step involves checking the quality of the downloaded sequencing data. The `fastqc` command is used to generate quality control reports for each of the fastq files. The `multiqc` command is then used to aggregate these reports into a single report.

```bash
fastqc *.fastq.gz
```
```
multiqc .
```

## Step 3: Create manifest file
This step involves creating a manifest file that lists all the fastq files along with their absolute paths. This file will be used in the next step for importing the data into QIIME 2.
```bash
ls -d "$PWD"/* > manifest.txt
```

## Step 4: Pack the paired-end data
This step involves importing the paired-end sequencing data into QIIME 2. The `qiime tools import` command is used for this purpose. The `--type 'SampleData[PairedEndSequencesWithQuality]'` option specifies the type of data being imported. The `--input-path manifest.csv` option specifies the path to the manifest file created in the previous step. The `--output-path demux` option specifies the output directory where the imported data will be stored. The `--input-format PairedEndFastqManifestPhred33` option specifies the format of the input data.
```bash
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest.csv --output-path demux --input-format PairedEndFastqManifestPhred33

```

## Step 5: Sequences quality visualisation (before the trimming)
This step involves visualizing the quality of the imported sequences before trimming. The `qiime demux summarize` command is used for this purpose. The `--i-data demux.qza` option specifies the input data, and the `--o-visualization demux` option specifies the output visualization.
```bash
qiime demux summarize --i-data demux.qza --o-visualization demux
```

## Step 6: Denoise Paired Sequences

This step involves denoising the paired-end sequences using DADA2, a pipeline for detecting and correcting (or removing) Illumina amplicon sequence data. The `qiime dada2 denoise-paired` command is used for this purpose.

This step is performed on 2023-06-02 using 430bp, 254-17 = 237, 234 - 21 = 213.

```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs primer-trimmed-demux_16s.qza \
  --p-trunc-len-f 327 \
  --p-trunc-len-r 213 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-max-ee-f 3 \
  --p-max-ee-r 3 \
  --p-n-threads 8 \
  --o-representative-sequences 430_327_213_rep-seqs.qza \
  --o-table 430_327_213_table.qza \
  --o-denoising-stats 430_327_213_stats.qza

```

## Step 7: Classify Sequences with Silva Classifier
This step involves classifying the sequences using a pre-trained Naive Bayes classifier and the Silva reference database. The `qiime feature-classifier classify-sklearn` command is used for this purpose.
```bash
qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/shared/scratch/pyau/qiime2-ref/silva-138-99-nb-classifier.qza \
  --i-reads 430_327_213_rep-seqs.qza \
  --p-n-jobs 4 \
  --o-classification 430_327_213_taxonomy.qza

```

## Step 8: Tabulate Metadata
This step involves tabulating metadata for easy visualization. The `qiime metadata tabulate` command is used for this purpose.
```bash
qiime metadata tabulate \
  --m-input-file 430_327_213_taxonomy.qza \
  --o-visualization 430_327_213_taxonomy.qzv

```

## Step 9: Taxonomy-based Filtering of Tables and Sequences
This step involves filtering tables and sequences based on taxonomy. Sequences that are classified as mitochondria or chloroplasts are excluded from further analysis as they do not represent microbial diversity. The `qiime taxa filter-table` and `qiime taxa filter-seqs` commands are used for this purpose.
```bash
qiime taxa filter-table \
  --i-table 430_327_213_table.qza \
  --i-taxonomy 430_327_213_taxonomy.qza \
  --p-include p__ \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table 430_327_213_table-with-phyla-no-mitochondria-no-chloroplast.qza
 ```
```
qiime taxa filter-seqs \
--i-sequences 430_327_213_rep-seqs.qza \
--i-taxonomy 430_327_213_taxonomy.qza \
--p-include p__ \
--p-exclude mitochondria,chloroplast \
--o-filtered-sequences 430_327_213_rep-seqs-with-phyla-no-mitochondria-no-chloroplast.qza
```