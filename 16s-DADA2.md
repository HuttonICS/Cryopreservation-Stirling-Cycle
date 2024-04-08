
# Bash Script for 16S Sequencing Data Analysis with Qiime2

This Bash script facilitates the analysis of 16S sequencing data using Qiime2. Ensure that necessary tools like fastqc, multiqc, and qiime are installed and properly configured in your environment. Also, replace the paths with actual paths where your files are located.

## Step 1: Download the 16S sequencing data from SRA database

Use the fastq-dump command from the SRA Toolkit to download sequencing data from the SRA database. This tool retrieves data in the FASTQ format, which is standard for storing biological sequences and their quality scores.

Ensure you have the 0SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit) installed. You can download it from the NCBI website, via Conda using the Bioconda channel, or through package managers like apt for Ubuntu.

```bash
fastq-dump --split-files --gzip
ERR6375874	ERR6375957	ERR6375997	ERR6376032
ERR6376063	ERR6376114	ERR6376145	ERR6376173
ERR6376240	ERR6376291	ERR6376334	ERR6376495
ERR6376576	ERR6376581	ERR6376586	ERR6376591
ERR6376592	ERR6376812	ERR6376813	ERR6376814
ERR6376815	ERR6377043	ERR6377045	ERR6377046
ERR6377048	ERR6377049	ERR6377050	ERR6377051
ERR6377052	ERR6377054	ERR6377055	ERR6377056
ERR6377058	ERR6377059	ERR6377060	ERR6377062
ERR6377137	ERR6377138	ERR6377139	ERR6377141
ERR6377142	ERR6377143	ERR6377144	ERR6377145
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

First, we will generate a `manifest.csv` file, which is an essential component in the process of creating a QIIME 2 artifact (.QZA) file. The `manifest.csv` file provides QIIME 2 with the necessary information about your sequencing reads. We also need to create a `manifest.txt` file that lists all directories in the current working directory. This can be done using the following command:

```bash
ls -d "$PWD"/* > manifest.txt
```

Next, we will use a script to transform the `manifest.txt` file into the `manifest.csv` format. You can copy the script below and paste it into a `.sh` file:

```bash
#!/bin/bash
# Input and output file paths
input_file="manifest.txt"
output_file="manifest.csv"

# Create the CSV header
echo "sample-id,absolute-filepath,direction" > "$output_file"

# Read the input file line by line
while IFS= read -r line; do
    # Extract the sample name, direction, and remove "_S123"
    sample=$(echo "$line" | sed 's#.*/##; s/\..*//; s/_L001_.*//; s/_\([0-9]\+\)$//; s/_S[0-9]\+$//')
    direction=$(echo "$line" | grep -o "_R[12]_")

    # Determine the direction and create the CSV lines
    if [ "$direction" = "_R1_" ]; then
        echo "$sample,$line,forward" >> "$output_file"
    elif [ "$direction" = "_R2_" ]; then
        echo "$sample,$line,reverse" >> "$output_file"
    fi
done < "$input_file"

echo "CSV file created: $output_file"
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
This step involves visualising the quality of the imported sequences before trimming. The `qiime demux summarize` command is used for this purpose. The `--i-data demux.qza` option specifies the input data, and the `--o-visualization demux` option specifies the output visualisation.

```bash
qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux
```

## Step 6: Denoise Paired Sequences

This step involves denoising the paired-end sequences using DADA2, a pipeline for detecting and correcting (or removing) Illumina amplicon sequence data. The `qiime dada2 denoise-paired` command is used for this purpose.

The number 430 represents the length of the amplified region targeting the V3-V4 hypervariable regions of the 16S rRNA gene. 

The parameter was estimated using FIGARO and can be found from https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/tree/main/16s/figaro

The calculations are as follows that the lenth of the primers need to be removed:

254-17 = 237

234 - 21 = 213

The parameters `--p-trunc-len-f` and `--p-trunc-len-r` for truncating sequences to a fixed length, ensuring high-quality data for downstream analysis. These parameters determine the maximum length to which forward and reverse reads will be truncated after quality filtering. Here, `--p-trunc-len-f 327` and `--p-trunc-len-r 213` indicate that forward reads will be truncated to 327 bases and reverse reads to 213 bases.

```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs primer-trimmed-demux_16s.qza \
  --p-trunc-len-f 237 \
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
