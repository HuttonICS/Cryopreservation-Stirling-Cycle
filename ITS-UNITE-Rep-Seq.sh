#!/bin/bash

### (1) Download the UNITE v8.3 Qiime2 formatted reference sequences and taxonomy files from https://unite.ut.ee/repository.php

### (2) Decompress the file using [tar xzf]  

### (3) Fix formatting errors
awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' sh_refs_qiime_ver8_dynamic_all_10.05.2021_dev.fasta | tr -d ' ' > sh_refs_qiime_ver8_dynamic_all_10.05.2021_dev_uppercase.fasta

### (4) Import the UNITE reference sequences into Qiime2
qiime tools import \
 --type FeatureData[Sequence] \
 --input-path sh_refs_qiime_ver8_dynamic_all_10.05.2021_dev_uppercase.fasta \
 --output-path unite-ver8-dynamic-seqs-10.05.2021.qza
 
### (5) Import the taxomony file
qiime tools import \
 --type FeatureData[Taxonomy] \
 --input-path sh_taxonomy_qiime_ver8_dynamic_all_10.05.2021_dev.txt \
 --output-path unite-ver8-dynamic-tax-10.05.2021.qza \
 --input-format HeaderlessTSVTaxonomyFormat

### (6) Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
 --i-reference-reads unite-ver8-dynamic-seqs-10.05.2021.qza \
 --i-reference-taxonomy unite-ver8-dynamic-tax-10.05.2021.qza \
 --o-classifier unite-ver8-dynamic-classifier-10.05.2021.qza
