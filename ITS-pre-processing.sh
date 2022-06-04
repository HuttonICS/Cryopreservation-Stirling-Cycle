#!/bin/bash

## (1) Download the ITS sequencing data from SRA database ##
fastq-dump --split-files --gzip ERR6377054	ERR6377055	ERR6377056	ERR6377058	ERR6377059	ERR6377060	ERR6377062	ERR6377137	ERR6377138	ERR6377139	ERR6377141	ERR6377142	ERR6377143	ERR6377144	ERR6377145

## The process used qiime2 Version 2022.2 on 28th May 2022 ##
## (2) applied fastqc and multiqc for sequence quality checking ##
fastqc *.gz
multiqc .

## (3) generate a manifest.txt and formatting for qiime2 ##
ls -d "$PWD"/* > manifest.txt

## (4) pack the paired-end data ##
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest.csv --output-path demux --input-format PairedEndFastqManifestPhred33

## (5) Sequences quality visualisation (before the trimming) ##
qiime demux summarize --i-data demux.qza --o-visualization demux

## (6) Cutadapt for the adaptors and primers ##
qiime cutadapt trim-paired --i-demultiplexed-sequences demux.qza --p-front-f TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGAACCWGCGGARGGATCA --p-front-r GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGCTGCGTTCTTCATCGATGC --p-error-rate 0 --o-trimmed-sequences trimmed-demux.qza --verbose > primer_trimming.log

## (7) Sequences quality visualisation (after the trimming) ##
qiime demux summarize --i-data trimmed-demux.qza --o-visualization trimmed-demux