#!/bin/bash

## denoise and merge using DADA2 ##
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimmed-demux.qza \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-n-threads 8 \
  --o-representative-sequences rep-seqs-unite.qza \
  --o-table table-unite.qza \
  --o-denoising-stats stats-unite.qza

## Map with UNITE reference database for Fungi taxonomy ## 
qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/shared/scratch/pyau/qiime2-ref/unite-8-dynamic-classifier-10.05.2021.qza \
  --i-reads rep-seqs-unite.qza \
  --o-classification taxonomy.qza

## Visualisation (Optional) ##
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
  
  ## Taxonomy-based filtering of tables and sequences for the unassigned sequences ##
  qiime taxa filter-table \
  --i-table table-unite.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude  Unassigned \
  --o-filtered-table table-unite-filtered.qza
  
  qiime taxa filter-seqs \
    --i-sequences rep-seqs-unite.qza \
    --i-taxonomy taxonomy.qza \
    --p-exclude Unassigned \
    --o-filtered-sequences rep-seqs-unite.filtered.qza
