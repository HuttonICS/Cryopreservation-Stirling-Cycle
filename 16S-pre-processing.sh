#!/bin/bash

## (1) Download the 16S sequencing data from SRA database ##
fastq-dump --split-files --gzip ERR6375874	ERR6375957	ERR6375997	ERR6376032	ERR6376063	ERR6376114	ERR6376145	ERR6376173	ERR6376240	ERR6376291	ERR6376334	ERR6376495	ERR6376576	ERR6376581	ERR6376586	ERR6376591	ERR6376592	ERR6376812	ERR6376813	ERR6376814	ERR6376815	ERR6377043	ERR6377045	ERR6377046	ERR6377048	ERR6377049	ERR6377050	ERR6377051	ERR6377052	ERR6377054	ERR6377055	ERR6377056	ERR6377058	ERR6377059	ERR6377060	ERR6377062	ERR6377137	ERR6377138	ERR6377139	ERR6377141	ERR6377142	ERR6377143	ERR6377144	ERR6377145

## (2) Quality Check ##
fastqc *.fastq.gz ## require fastqc
multiqc . ## require multiqc

## (3) Create manifest file ##
ls -d "$PWD"/* > manifest.txt

## (4) pack the paired-end data ##
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest.csv --output-path demux --input-format PairedEndFastqManifestPhred33

## (5) Sequences quality visualisation (before the trimming) ##
qiime demux summarize --i-data demux.qza --o-visualization demux